import argparse
import traceback
import os
import errno
import pathlib
import logging
import numpy
import sys

"""
Outputs all genes that passed a filter and their associated counts.
Takes two files as input.
Input file 1: the names of the genes that passed a filter.
Input file 2: genes with their counts before applying filter.

TO DO: Probably no need to use numpy here. Simpler with dict?
TO DO: Replace Java style for loops with Python style.
TO DO: Automatically detect if input file is tsv, csv, ssv.
"""
class FilterCounts(object):

    global dictvalue
    #Variable used to show a key exists in the dict to replace None
    dictvalue = 1

    def __init__(self, input_filtered_genes, input_genes_counts, output, debug):
        filename = os.path.basename(__file__)
        f = filename.split(".")
        prefix_filename = f[0]
        self.input_filtered_genes = input_filtered_genes
        self.input_genes_counts = input_genes_counts
        self.debug = debug
        self.log_filename = prefix_filename + '.log'
        self.logger = self.get_logger(self.log_filename)
        self.logger.info('----------------------------------------------')
        self.logger.info('Running: ' + filename)
        self.logger.info('Input file parameter 1: ' + self.input_filtered_genes)
        self.logger.info('Input file parameter 2: ' + self.input_genes_counts)
        self.output = output
        self.DELIMITER=','  # TO DO: make this settable

    def assert_file_exists (self, filename):
        if not os.path.isfile(filename):
            self.logger.error('File not found: '+filename)
            raise FileNotFoundError(filename)

    def __del__(self):
        logging.shutdown()
        if self.logger and self.logger.handlers:
            self.logger.handlers.clear()

    def filter_counts(self):
        dict_filtered = self.file_to_dict(self.input_filtered_genes)
        array_counts = self.file_to_array(self.input_genes_counts)
        array_counts_filtered = self.apply_filter(dict_filtered, array_counts)
        if (self.output != None):
            numpy.savetxt(self.output, array_counts_filtered, fmt='%s')
            self.logger.info('Output file: ' + self.output)
            self.assert_file_exists(self.output)
        else:
            numpy.savetxt(sys.stdout, array_counts_filtered, fmt='%s')

    def file_to_array(self, input):
        """Converts a file to a 1D array with each line of the file a different value in the array"""
        self.assert_file_exists(input)
        a = []
        infile = open(input, 'r')
        for line in infile:
            line = line.strip() #removes \n for new line from the values
            a.append(line)
        infile.close()
        return a

    def apply_filter(self, dict_filtered, array_counts):
        filter_array = []
        for r in range(len(array_counts)):
            line = array_counts[r]
            key = line.split(self.DELIMITER, 1) [0]
            if key in dict_filtered:
                filter_array.append(line)
        return(filter_array)

    def file_to_dict(self, filename):
        self.assert_file_exists(filename)
        dict = {}
        infile = open(filename, 'r')
        for line in infile:
            key = line.strip()
            d = {key : dictvalue}
            dict.update(d)
        infile.close()
        return dict

    def get_logger(self, filename):
        logger = logging.getLogger(filename)
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter(
            '%(asctime)s: %(levelname)s: %(message)s')
        file_handler = logging.FileHandler(filename)
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        return logger

    def main(self):
        self.filter_counts()

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description = 'Filters the collation of three replicates.')
    parser.add_argument('input_filtered_genes', help='File of genes that passed  filter',
        type=str)
    parser.add_argument('input_genes_counts', help='File containing genes with their counts (*.csv)',
        type=str)
    parser.add_argument('--output', help='Name for output file that contains filtered genes with their counts',
        type=str)
    parser.add_argument('--debug', action='store_true', help='print traceback')
    args = parser.parse_args()

if __name__ =="__main__":
    try:
        args_parse()
        instance = FilterCounts(args.input_filtered_genes, args.input_genes_counts, args.output, args.debug)
        instance.main()
    except Exception as e:
        instance.logger.error(e)
        print('ERROR! check ' + instance.log_filename)
        if args.debug:
            print(traceback.format_exc())
        else:
            print('ERROR, run with --debug for traceback')
        exit(1)
