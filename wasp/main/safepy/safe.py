#! /usr/bin/env python
"""This file contains the code for the SAFE class and command-line access."""

import configparser
import os
from pathlib import Path
import sys
import textwrap
import argparse
import pickle
import time
import re
import logging

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import multiprocessing as mp

from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import hypergeom
from itertools import compress
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import fdrcorrection

from .safe_io import *
from .safe_extras import *

#from . import safe_io
#from . import safe_extras


class SAFE:
    """
    Defines an instance of SAFE analysis.
    Contains all data, all parameters and provides the main methods for performing analysis.
    """

    def __init__(self,
                 path_to_ini_file='',
                 path_to_safe_data=None,
                 verbose=True):
        """
        Initiate a SAFE instance and define the main settings for analysis.
        The settings are automatically extracted from the specified (or default) INI configuration file.
        Alternatively, each setting can be changed manually after initiation.

        :param path_to_ini_file (str): Path to the configuration file. If not specified, safe_default.ini will be used.
        :param verbose (bool): Defines whether or not intermediate output will be printed out.

        """

        self.verbose = verbose

        self.default_config = None

        self.path_to_safe_data = path_to_safe_data
        self.path_to_network_file = None
        self.view_name = None
        self.path_to_attribute_file = None

        self.graph = None
        self.node_key_attribute = 'label_orf'

        self.attributes = None
        self.nodes = None
        self.node2attribute = None
        self.num_nodes_per_attribute = None
        self.attribute_sign = 'both'

        self.node_distance_metric = 'shortpath_weighted_layout'
        self.neighborhood_radius_type = None
        self.neighborhood_radius = None

        self.background = 'attribute_file'
        self.num_permutations = 1000
        self.multiple_testing = False
        self.neighborhood_score_type = 'sum'
        self.enrichment_type = 'auto'
        self.enrichment_threshold = 0.05
        self.enrichment_max_log10 = 16
        self.attribute_enrichment_min_size = 10

        self.neighborhoods = None

        self.ns = None
        self.pvalues_neg = None
        self.pvalues_pos = None
        self.nes = None
        self.nes_threshold = None
        self.nes_binary = None

        self.attribute_unimodality_metric = 'connectivity'
        self.attribute_distance_metric = 'jaccard'
        self.attribute_distance_threshold = 0.75

        self.domains = None
        self.node2domain = None

        # Output
        self.output_dir = ''

        # Read both default and user-defined settings
        self.read_config(path_to_ini_file, path_to_safe_data=self.path_to_safe_data)

        # Validate config
        self.validate_config()

    def read_config(self, path_to_ini_file, path_to_safe_data=None):

        """
        Read the settings from an INI file and update the attributes in the SAFE class.

        :param path_to_ini_file (str): Path to the configuration file. If not specified, safe_default.ini will be used.
        :return: none
        """

        # Location of this code
        loc = os.path.dirname(os.path.abspath(__file__))

        # Load default settings
        default_config_path = os.path.join(loc, 'safe_default.ini')
        default_config = configparser.ConfigParser(allow_no_value=True,
                                                   comment_prefixes=('#', ';', '{'),
                                                   inline_comment_prefixes='#')

        with open(default_config_path, 'r') as f:
            default_config.read_file(f)

        self.default_config = default_config['DEFAULT']

        # Load user-defined settings, if any
        config = configparser.ConfigParser(defaults=default_config['DEFAULT'],
                                           allow_no_value=True,
                                           comment_prefixes=('#', ';', '{'),
                                           inline_comment_prefixes='#')
        config.read(path_to_ini_file)

        if 'Input files' not in config:
            config['Input files'] = {}
        if path_to_safe_data is None:
            path_to_safe_data = config.get('Input files', 'safe_data')  # falls back on default if empty
            if path_to_safe_data == '':
                path_to_safe_data = None
        path_to_network_file = config.get('Input files', 'networkfile')  # falls back on default if empty
        path_to_attribute_file = config.get('Input files', 'annotationfile')  # falls back on default if empty

        self.path_to_safe_data = path_to_safe_data
        if not self.path_to_safe_data is None:
            assert self.path_to_safe_data.endswith('/'), "path_to_safe_data should end with '/', else `os.path.join` may not provide desired output."
            self.path_to_network_file = os.path.join(self.path_to_safe_data, path_to_network_file)
            self.path_to_attribute_file = os.path.join(self.path_to_safe_data, path_to_attribute_file)
        else:
            ## direct paths to the network and attribute files
            self.path_to_network_file = path_to_network_file
            self.path_to_attribute_file = path_to_attribute_file
            
        self.attribute_sign = config.get('Input files', 'annotationsign') # falls back on default if empty

        if 'Analysis parameters' not in config:
            config['Analysis parameters'] = {}

        self.background = config.get('Analysis parameters', 'background')
        self.node_distance_metric = config.get('Analysis parameters', 'nodeDistanceType')
        self.neighborhood_radius_type = config.get('Analysis parameters', 'neighborhoodRadiusType')
        self.neighborhood_radius = float(config.get('Analysis parameters', 'neighborhoodRadius'))

        self.attribute_unimodality_metric = config.get('Analysis parameters', 'unimodalityType')
        self.attribute_distance_metric = config.get('Analysis parameters', 'groupDistanceType')
        self.attribute_distance_threshold = float(config.get('Analysis parameters', 'groupDistanceThreshold'))

        self.output_dir = os.path.dirname(path_to_ini_file)
        if not self.output_dir:
            self.output_dir = loc

    def validate_config(self):

        """
        Test the validity of the current settings in the SAFE class before running the analysis.

        :return: none
        """

        # Check that the option parameters are valid
        if self.background not in ['attribute_file', 'network']:
            user_setting = self.background
            self.background = self.default_config.get('background')    # Restore the default value.
            raise ValueError(('%s is not a valid setting for background. '
                              'Valid options are: attribute_file, network.' % user_setting))

        if self.node_distance_metric not in ['euclidean', 'shortpath', 'shortpath_weighted_layout']:
            user_setting = self.node_distance_metric
            self.node_distance_metric = self.default_config.get('nodeDistanceType')    # Restore the default value.
            raise ValueError(('%s is not a valid setting for node_distance_metric. '
                              'Valid options are: euclidean, shortpath, shortpath_weighted_layout' % user_setting))

        if self.attribute_sign not in ['highest', 'lowest', 'both']:
            user_setting = self.attribute_sign
            self.attribute_sign = self.default_config.get('annotationsign')   # Restore the default value.
            raise ValueError(('%s is not a valid setting for attribute_sign. '
                              'Valid options are: highest, lowest, both' % user_setting))

        if not isinstance(self.num_permutations, int) or (self.num_permutations < 10):
            self.num_permutations = 1000    # Restore the default value.
            raise ValueError('num_permutations must be an integer equal or greater than 10.')

        if not isinstance(self.enrichment_threshold, float) or (self.enrichment_threshold <= 0) or (self.enrichment_threshold >= 1):
            self.enrichment_threshold = 0.05    # Restore the default value.
            raise ValueError('enrichment_threshold must be in the (0,1) range.')

        if not isinstance(self.enrichment_max_log10, (int, float)):
            self.enrichment_max_log10 = 16    # Restore the default value.
            raise ValueError('enrichment_max_log10 must be a number.')

        if not isinstance(self.attribute_enrichment_min_size, int) or (self.attribute_enrichment_min_size < 2):
            self.attribute_enrichment_min_size = 10    # Restore the default value.
            raise ValueError('attribute_enrichment_min_size must be an integer equal or greater than 2.')

        if not isinstance(self.attribute_distance_threshold, float) or (self.attribute_distance_threshold <= 0) or (self.attribute_distance_threshold >= 1):
            self.attribute_distance_threshold = 0.75    # Restore the default value.
            raise ValueError('attribute_enrichment_min_size must be a float number in the (0,1) range.')


    def load_network(self, **kwargs):
        """
        Load the network from a source file and, if necessary, apply a network layout.

        Keyword Args:
            * network_file (:obj:`str`, optional): Path to the file containing the network. Note: if the path to safe data (`path_to_safe_data`) is provided, this would the path inside the `safe_data` folder, else a direct path to the file. 
            * node_key_attribute (:obj:`str`, optional): Name of the node attribute that should be treated as key identifier.

        :return: none
        """

        # Overwriting the global settings, if required
        if 'network_file' in kwargs:
            if self.path_to_safe_data is None:
                self.path_to_network_file = kwargs['network_file']
            else:
                self.path_to_network_file = os.path.join(self.path_to_safe_data, kwargs['network_file'])
            del kwargs['network_file'] ## remove the redundant/old path
        assert os.path.exists(self.path_to_network_file), self.path_to_network_file # os.path.join may misbehave if there are extra '/' at the place where the paths are joined.
        if 'view_name' in kwargs:
            self.view_name = kwargs['view_name']
        if 'node_key_attribute' in kwargs:
            self.node_key_attribute = kwargs['node_key_attribute']

        # Make sure that the settings are still valid
        self.validate_config()

        if type(self.path_to_network_file) == nx.Graph:

            self.graph = self.path_to_network_file

        else:

            # [_, file_extension] = os.path.splitext(self.path_to_network_file)
            file_extension=Path(self.path_to_network_file).suffixes[0] # compatible with double extension e.g. txt.gz
            if self.verbose:
                logging.info('Loading network from %s' % self.path_to_network_file)

            if file_extension in ['.txt','.tsv']:
                self.graph = load_network_from_txt(self.path_to_network_file,
                                                   node_key_attribute=self.node_key_attribute,
                                                   verbose=self.verbose)

        # Setting the node key for mapping attributes
        key_list = nx.get_node_attributes(self.graph, self.node_key_attribute)

        if not bool(key_list):
            raise Exception('The specified node key attribute (%s) does not exist in this network. '
                            'These attributes exist instead: %s. '
                            'Set node_key_attribute to one of these options.'
                            % (self.node_key_attribute, ', '.join(self.graph.node[0].keys())))
        else:
            nx.set_node_attributes(self.graph, key_list, name='key')
            label_list = nx.get_node_attributes(self.graph, 'label')
            self.nodes = pd.DataFrame(data={'id': list(label_list.keys()),
                                            'key': list(key_list.values()),
                                            'label': list(label_list.values())})


    def load_attributes(self, **kwargs):
        """
        Preprocess and load the attributes i.e. features of the genes.
        
        Keyword arguments:
            kwargs: parameters provided to `read_attributes` function.
            * attribute_file (:obj:`str`, optional): Path to the file containing the attributes. Note: if path to safe data (`path_to_safe_data`) is provided, this would the path inside the `safe_data` folder, else a direct path to the file. 
        """
        
        # Overwrite the global settings, if required
        if 'attribute_file' in kwargs:
            if self.path_to_safe_data is None or isinstance(kwargs['attribute_file'],pd.DataFrame):
                self.path_to_attribute_file = kwargs['attribute_file']
            elif isinstance(kwargs['attribute_file'],str):
                self.path_to_attribute_file = os.path.join(self.path_to_safe_data, kwargs['attribute_file'])
            else:
                raise ValueError(type(kwargs['attribute_file']))     
            del kwargs['attribute_file'] ## remove the redundant/old path
        if isinstance(self.path_to_attribute_file,str):
            assert os.path.exists(self.path_to_attribute_file), self.path_to_attribute_file # os.path.join may misbehave if there are extra '/' at the place where the paths are joined.
            
        # Make sure that the settings are still valid
        self.validate_config()

        node_label_order = list(nx.get_node_attributes(self.graph, self.node_key_attribute).values())

        if self.verbose and isinstance(self.path_to_attribute_file, str):
            logging.info('Loading attributes from %s' % self.path_to_attribute_file)

        [self.attributes, _, self.node2attribute] = read_attributes(node_label_order=node_label_order,
                                                                    verbose=self.verbose, 
                                                                    attribute_file=self.path_to_attribute_file, 
                                                                    **kwargs)

    def define_neighborhoods(self, **kwargs):
        """
        
        """
        # Overwriting the global settings, if required
        if 'node_distance_metric' in kwargs:
            self.node_distance_metric = kwargs['node_distance_metric']

        if 'neighborhood_radius_type' in kwargs:
            self.neighborhood_radius_type = kwargs['neighborhood_radius_type']

        if 'neighborhood_radius' in kwargs:
            self.neighborhood_radius = kwargs['neighborhood_radius']

        # Make sure that the settings are still valid
        self.validate_config()

        all_shortest_paths = {}
        neighborhoods = np.zeros([self.graph.number_of_nodes(), self.graph.number_of_nodes()], dtype=int)

        if self.node_distance_metric == 'euclidean':
            x = list(dict(self.graph.nodes.data('x')).values())
            nr = self.neighborhood_radius * (np.max(x) - np.min(x))

            x = np.matrix(self.graph.nodes.data('x'))[:, 1]
            y = np.matrix(self.graph.nodes.data('y'))[:, 1]

            node_coordinates = np.concatenate([x, y], axis=1)
            node_distances = squareform(pdist(node_coordinates, 'euclidean'))

            neighborhoods[node_distances < nr] = 1

        else:

            if self.node_distance_metric == 'shortpath_weighted_layout':
                # x = np.matrix(self.graph.nodes.data('x'))[:, 1]
                x = list(dict(self.graph.nodes.data('x')).values())

                
                nr = self.neighborhood_radius * (np.max(x) - np.min(x))
                all_shortest_paths = dict(nx.all_pairs_dijkstra_path_length(self.graph,
                                                                            weight='length', cutoff=nr))
            elif self.node_distance_metric == 'shortpath':
                nr = self.neighborhood_radius
                all_shortest_paths = dict(nx.all_pairs_dijkstra_path_length(self.graph, cutoff=nr))

            neighbors = [(s, t) for s in all_shortest_paths for t in all_shortest_paths[s].keys()]

            for i in neighbors:
                neighborhoods[i] = 1

            self.node_distances = all_shortest_paths

        # Set diagonal to zero (a node is not part of its own neighborhood)
        # np.fill_diagonal(neighborhoods, 0)

        # Calculate the average neighborhood size
        num_neighbors = np.sum(neighborhoods, axis=1)

        if self.verbose:
            logging.info('Node distance metric: %s' % self.node_distance_metric)
            logging.info('Neighborhood definition: %.2f x %s' % (self.neighborhood_radius, self.neighborhood_radius_type))
            logging.info('Number of nodes per neighborhood (mean +/- std): %.2f +/- %.2f' % (np.mean(num_neighbors), np.std(num_neighbors)))

        self.neighborhoods = neighborhoods

    def compute_pvalues(self, **kwargs):

        if 'how' in kwargs:
            self.enrichment_type = kwargs['how']

        if 'neighborhood_score_type' in kwargs:
            self.neighborhood_score_type = kwargs['neighborhood_score_type']

        if 'multiple_testing' in kwargs:
            self.multiple_testing = kwargs['multiple_testing']

        if 'background' in kwargs:
            self.background = kwargs['background']

        # Make sure that the settings are still valid
        self.validate_config()

        if self.background == 'network':
            logging.info('Setting all null attribute values to 0. Using the network as background for enrichment.')
            self.node2attribute[np.isnan(self.node2attribute)] = 0

        num_vals = self.node2attribute.shape[0]
        num_nans = np.sum(np.isnan(self.node2attribute), axis=0)

        if any(num_nans/num_vals > 0.5):
            logging.warning("WARNING: more than 50% of nodes in the network as set to NaN and will be ignored for calculating enrichment.\n'Consider setting sf.background = ''network''.'")

        # Warn users if more than 50% of values are NaN
        num_other_values = np.sum(~np.isnan(self.node2attribute) & ~np.isin(self.node2attribute, [0, 1]))

        if (self.enrichment_type == 'hypergeometric') or ((self.enrichment_type == 'auto') and (num_other_values == 0)):
            self.compute_pvalues_by_hypergeom(**kwargs)
        else:
            self.compute_pvalues_by_randomization(**kwargs)

        idx = ~np.isnan(self.nes)
        self.nes_binary = np.zeros(self.nes.shape)
        self.nes_binary[idx] = np.abs(self.nes[idx]) > -np.log10(self.enrichment_threshold)

        self.attributes['num_neighborhoods_enriched'] = np.sum(self.nes_binary, axis=0)

    def compute_pvalues_by_randomization(self, **kwargs):

        if kwargs:
            logging.warning('Current settings (possibly overwriting global ones):')
            for k in kwargs:
                logging.warning('\t%s=%s' % (k, str(kwargs[k])))

        logging.info('Using randomization to calculate enrichment...')

        # Pause for 1 sec to prevent the progress bar from showing up too early
        time.sleep(1)

        if 'num_permutations' in kwargs:
            self.num_permutations = kwargs['num_permutations']

        num_processes = 1
        if 'processes' in kwargs:
            num_processes = kwargs['processes']

        # Make sure that the settings are still valid
        self.validate_config()

        N_in_neighborhood_in_group = compute_neighborhood_score(self.neighborhoods,
                                                                self.node2attribute,
                                                                self.neighborhood_score_type)
        self.ns = N_in_neighborhood_in_group

        if num_processes > 1:

            num_permutations_x_process = np.ceil(self.num_permutations / num_processes).astype(int)
            self.num_permutations = num_permutations_x_process * num_processes

            arg_tuple = (self.neighborhoods, self.node2attribute,
                         self.neighborhood_score_type, num_permutations_x_process)
            list_for_parallelization = [arg_tuple] * num_processes

            ctx = mp.get_context('spawn')
            pl = ctx.Pool(processes=num_processes)
            res = pl.map(run_permutations, list_for_parallelization)
            pl.close()
            pl.join()

            [counts_neg_list, counts_pos_list] = map(list, zip(*res))

            counts_neg = np.sum(np.stack(counts_neg_list, axis=2), axis=2)
            counts_pos = np.sum(np.stack(counts_pos_list, axis=2), axis=2)

        else:

            arg_tuple = (self.neighborhoods, self.node2attribute,
                         self.neighborhood_score_type, self.num_permutations)
            [counts_neg, counts_pos] = run_permutations(arg_tuple)

        idx = np.isnan(N_in_neighborhood_in_group)
        counts_neg[idx] = np.nan
        counts_pos[idx] = np.nan

        self.pvalues_neg = counts_neg / self.num_permutations
        self.pvalues_pos = counts_pos / self.num_permutations

        # Correct for multiple testing
        if self.multiple_testing:
            logging.info('Running FDR-adjustment of p-values...')
            out = np.apply_along_axis(fdrcorrection, 1, self.pvalues_neg)
            self.pvalues_neg = out[:, 1, :]

            out = np.apply_along_axis(fdrcorrection, 1, self.pvalues_pos)
            self.pvalues_pos = out[:, 1, :]

        # Log-transform into neighborhood enrichment scores (NES)
        # Necessary conservative adjustment: when p-value = 0, set it to 1/num_permutations
        nes_pos = -np.log10(np.where(self.pvalues_pos == 0, 1/self.num_permutations, self.pvalues_pos))
        nes_neg = -np.log10(np.where(self.pvalues_neg == 0, 1/self.num_permutations, self.pvalues_neg))

        if self.attribute_sign == 'highest':
            self.nes = nes_pos
        elif self.attribute_sign == 'lowest':
            self.nes = nes_neg
        elif self.attribute_sign == 'both':
            self.nes = nes_pos - nes_neg

    def compute_pvalues_by_hypergeom(self, **kwargs):

        if kwargs:
            if 'verbose' in kwargs:
                self.verbose = kwargs['verbose']

            if self.verbose:
                logging.warning('Overwriting global settings:')
                for k in kwargs:
                    logging.warning('\t%s=%s' % (k, str(kwargs[k])))

        # Make sure that the settings are still valid
        self.validate_config()

        if self.verbose:
            logging.info('Using the hypergeometric test to calculate enrichment...')

        # Nodes with not-NaN values in >= 1 attribute
        nodes_not_nan = np.any(~np.isnan(self.node2attribute), axis=1)

        # -- Number of nodes
        # n = self.graph.number_of_nodes()    # total
        n = np.sum(nodes_not_nan)    # with not-NaN values in >=1 attribute

        N = np.zeros([self.graph.number_of_nodes(), len(self.attributes)]) + n

        # -- Number of nodes annotated to each attribute
        N_in_group = np.tile(np.nansum(self.node2attribute, axis=0), (self.graph.number_of_nodes(), 1))

        # -- Number of nodes in each neighborhood
        # neighborhood_size = np.sum(self.neighborhoods, axis=0)[:, np.newaxis]    # total
        neighborhood_size = np.dot(self.neighborhoods,
                                   nodes_not_nan.astype(int))[:, np.newaxis] # with not-NaN values in >=1 attribute

        N_in_neighborhood = np.tile(neighborhood_size, (1, len(self.attributes)))

        # -- Number of nodes in each neighborhood and  annotated to each attribute
        N_in_neighborhood_in_group = np.dot(self.neighborhoods,
                                            np.where(~np.isnan(self.node2attribute), self.node2attribute, 0))

        self.pvalues_pos = hypergeom.sf(N_in_neighborhood_in_group - 1, N, N_in_group, N_in_neighborhood)

        # Correct for multiple testing
        if self.multiple_testing:

            if self.verbose:
                logging.info('Running FDR-adjustment of p-values...')

            out = np.apply_along_axis(fdrcorrection, 1, self.pvalues_pos)
            self.pvalues_pos = out[:, 1, :]

        # Log-transform into neighborhood enrichment scores (NES)
        self.nes = -np.log10(self.pvalues_pos)


    def plot_network(self, 
                     foreground_color = '#ffffff',
                     background_color='#000000',
                     labels=[],
                     **kwargs_mark_nodes,
                    ):
        """
        Plot the base network.
        
        Parameters:
            labels (list): the genes to show on the network.
        
        Keyword parameters:
            kwargs_mark_nodes: parameters provided to `mark_nodes` function.
        """
        ax=plot_network(self.graph, background_color=background_color)
        # Plot the labels, if any
        if len(labels)!=0:
            ## get the coordinates of the points
            node_xy_labels,labels_found=get_node_coordinates(graph=self.graph,labels=labels)
            ## mark the nodes
            ax=mark_nodes(
                       x=node_xy_labels[:, 0],
                       y=node_xy_labels[:, 1],
                       labels=labels_found,
                       ax=ax,
                       foreground_color=foreground_color,
                       background_color=background_color,
                       **kwargs_mark_nodes,
                      )
        return ax
    

    def print_output_files(self, **kwargs):

        if 'output_dir' in kwargs:
            self.output_dir = kwargs['output_dir']

        # Domain properties
        path_domains = f"{self.output_dir}domain_properties_annotation.txt"
        if self.domains is not None:
            self.domains.drop(labels=[0], axis=0, inplace=True, errors='ignore')
            self.domains.to_csv(path_domains, sep='\t')
            logging.info(path_domains)

        # Attribute properties
        path_attributes = f"{self.output_dir}attribute_properties_annotation.txt"
        self.attributes.to_csv(path_attributes, sep='\t')
        logging.info(path_attributes)

        # Node properties
        path_nodes = f"{self.output_dir}node_properties_annotation.txt"

        t = nx.get_node_attributes(self.graph, 'key')
        ids = list(t.keys())
        keys = list(t.values())
        t = nx.get_node_attributes(self.graph, 'label')
        labels = list(t.values())
        if self.node2domain is not None:
            domains = self.node2domain['primary_domain'].values
            ness = self.node2domain['primary_nes'].values
            num_domains = self.node2domain[self.domains['id']].sum(axis=1).values
            self.nodes = pd.DataFrame(data={'id': ids, 'key': keys, 'label': labels, 'domain': domains,
                                            'nes': ness, 'num_domains': num_domains})
        else:

            self.nodes = pd.DataFrame(self.nes)
            self.nodes.columns = self.attributes['name']
            self.nodes.insert(loc=0, column='key', value=keys)
            self.nodes.insert(loc=1, column='label', value=labels)

        self.nodes.to_csv(path_nodes, sep='\t')
        logging.info(path_nodes)


def run_safe_batch(attribute_file):

    sf = SAFE()
    sf.load_network()
    sf.define_neighborhoods()

    sf.load_attributes(attribute_file=attribute_file)
    sf.compute_pvalues(num_permutations=1000)

    return sf.nes


if __name__ == '__main__':

    start = time.time()

    parser = argparse.ArgumentParser(description='Run Spatial Analysis of Functional Enrichment (SAFE) on the default Costanzo et al., 2016 network')
    parser.add_argument('path_to_attribute_file', metavar='path_to_attribute_file', type=str,
                        help='Path to the file containing label-to-attribute annotations')

    args = parser.parse_args()

    # Load the attribute file
    [attributes, node_label_order, node2attribute] = read_attributes(args.path_to_attribute_file)

    nr_processes = mp.cpu_count()
    nr_attributes = attributes.shape[0]

    chunk_size = np.ceil(nr_attributes / nr_processes).astype(int)
    chunks = np.array_split(np.arange(nr_attributes), nr_processes)

    all_chunks = []
    for chunk in chunks:
        this_chunk = pd.DataFrame(data=node2attribute[:, chunk], index=node_label_order,
                                  columns=attributes['name'].values[chunk])
        all_chunks.append(this_chunk)

    pool = mp.Pool(processes=nr_processes)

    combined_nes = []

    logging.info('Running SAFE on %d chunks of size %d...' % (nr_processes, chunk_size))
    for res in pool.map_async(run_safe_batch, all_chunks).get():
        combined_nes.append(res)

    all_nes = np.concatenate(combined_nes, axis=1)

    output_file = format('%s_safe_nes.p' % args.path_to_attribute_file)

    logging.info('Saving the results...')
    with open(output_file, 'wb') as handle:
        pickle.dump(all_nes, handle)

