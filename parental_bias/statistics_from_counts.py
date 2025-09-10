# See MOLBAR technical specification TS3
"""
This module generates statistics from read counts.
This uses limma in R to generate t-stat, P-value, adjusted P-value.
Read counts are accepted in a space-delimited text file.
(To do: Accept csv also.)
Expected input has 3 replicates of (maternal,paternal) counts.
(To do: Accept any number >= 3.)
Statistics are given in a csv file.
This writes a log file (or appends to the existing one).

Environment setup:
$ R
> install.packages(“BiocManager”)
> library(BiocManager)
> install(“limma”)
> install(“gtools”)
> quit()
$
Environment setup on Saga:
(Unfortunately must use version "free open source software 2018b"
to match the only available version of RSEM on Saga.)
$ module load R/3.5.1-foss-2018b
$ module load  R-bundle-Bioconductor/3.8-foss-2018b-R-3.5.1
"""

import argparse
import os
import errno
# import magic
import logging
import json
# import re -unused library
# from operator import itemgetter -unused library
# import pydoc- unused library but may be useful later
import shutil
import datetime
#import subprocess

class StatisticsFromCounts():
    """
    A class used to run a dataset through an R script and then produce
    a final csv file from the script output and the dataset.
    This is a wrapper for the limma.foldchange.r script.

    Typical usage:
    >>> stats = statistics_from_counts.StatisticsFromCounts(
    ...         script_filename,input_filename,True)
    >>> stats.run_all()
    """

    # class variables
    logger = logging.getLogger()
    COMMAND = None
    R_SCRIPT_RUNNER = "Rscript"
    REMOVE_HEADERS = True

    def __init__(self, script_, cross_, debug_=False):
        """
        script      [<path>/]<filename> e.g. 'limma.foldchange.r'
        cross       [<path>/]<filename> e.g. 'data/Col_x_Ler.read_counts.txt'
        debug     Write progress to log (default = write warnings & errors)
        """
        self.script = script_
        self.cross = cross_
        self.BASENAME = os.path.basename(self.cross)
        self.sorted_stats_file = self.BASENAME + ".de.sorted"
        self.unsorted_filename = self.BASENAME + ".de"
        self.results_filename = self.BASENAME + ".final.csv"
        self.myrec = {''}
        self.raw_counts = ['']
        if debug_:
            loglevel = logging.INFO
        else:
            loglevel = logging.WARNING
        logfilename = 'statistics_from_counts.log'
        logging.basicConfig(
            filename=logfilename,
            level=loglevel,
            format='%(asctime)s: %(levelname)s: %(message)s')
        self.assert_file_exists('script file', self.script)
        self.assert_file_exists('cross file', self.cross)
        self.script_check(self.script)

    def assert_file_exists(self, descrip, filename):
        rlogger = StatisticsFromCounts.logger
        if not os.path.isfile(filename):
            rlogger.error(
                "%s and (%s) not found.", descrip, filename)
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), filename)

    def load_raw_counts(self):
        """Load the counts associated with each gene.
        Assume the data file uses a single-space delimiter.
        Get ready to write this in csv format later.
        Replace spaces with commas and add trailing comma.
        """
        with open(self.cross, 'r') as myfile:
            # Read entire file into a string.
            # Global replace space with comma.
            # Split the string into a list of strings.
            list1 = myfile.read().replace(' ', ',').split('\n')
            myfile.close()
        # Use python list comprehension to add comma to end of every line.
        self.raw_counts = [x + ',' for x in list1]
        # Fix this later: we always get an extra (blank) entry.
        # Downstream code breaks if we remove it.
        # For now, just remove comma from the extra entry.
        if self.raw_counts[-1] == ',':
            self.raw_counts[-1] = ''

    def merge_raw_stats(self):
        """Load the stats associated with each gene.
        Assume the data file uses a comma delimiter.
        Assume the number of lines matches raw counts.
        Remove double quotes.
        Discard the first line, assumed to hold column headers.
        Use first number as index into raw_counts.
        Merge raw stats with raw counts and write the file.
        """
        with open(self.unsorted_filename, 'r') as myfile:
            if StatisticsFromCounts.REMOVE_HEADERS:
                first_data_line = 1   # currently required
            else:
                first_data_line = 0   # in case things change
            index_off_by_one = 1
            line_num = 0
            for line1 in myfile:
                if line_num >= first_data_line:
                    # Global remove double quotes.
                    line2 = line1.replace('"', '')
                    first_comma = line2.find(',')
                    gene_num = 0
                    if first_comma > 0:
                        gene_num = int(line2[0:first_comma])
                    if gene_num > 0:
                        # Gene "1" was at position 0 in raw_counts
                        self.raw_counts[gene_num-index_off_by_one] += line2
                line_num += 1
            myfile.close()
        with open(self.results_filename, 'w') as outfile:
            for big_str in self.raw_counts:
                outfile.write(big_str)
            outfile.close()
        self.assert_file_exists('R result file', self.results_filename)

    def remove_file(self, filename):
        rlogger = StatisticsFromCounts.logger
        try:
            os.remove(filename)
        except OSError as e:
            rlogger.warning(filename+" not found and not removed. \
                            Error: %s - %s." % (e.filename, e.strerror))

    def script_check(self,filename):
        rlogger = StatisticsFromCounts.logger
        script_ext = os.path.splitext(filename)[1]
        if script_ext != ".r" and script_ext != ".R":
            rlogger.warning(
                "Script argument " + filename +
                "should have a .r extension. Yours has " +
                script_ext)

    def script_run(self):
        """Runs the R script on the dataset"""
        # Note the R script writes unsorted temp file <infile>.de
        rlogger = StatisticsFromCounts.logger
        if self.cross != self.BASENAME:
            # Copy input file to local directory
            # so the R script will read and write here.
            # This way, output files go to current directory.
            shutil.copy(self.cross, self.BASENAME)
        COMMAND = "{} {} {}".format(StatisticsFromCounts.R_SCRIPT_RUNNER,
                                    self.script,
                                    self.BASENAME)
        rlogger.info("Compute statistics on %s", self.BASENAME)
        rlogger.info("Getting a result file requires the 'Rscript' program on \
                      your system.")
        rlogger.info(COMMAND)
        # The following command may raise an exception.
        # However, it may exit cleanly with no output file.
        # TO DO: help user figure out what went wrong here.
        #subprocess.call(
        #    [StatisticsFromCounts.R_SCRIPT_RUNNER,
        #     self.script,
        #     self.BASENAME])
        os.system(COMMAND)
        self.assert_file_exists('R result file',
                                self.unsorted_filename)

    def run_recorder(self):
        """Records information about the run to a file"""
        module_name = __file__
        date_str = str(datetime.datetime.now())
        self.myrec = {'prog-name': module_name,
                      'input_files': [self.cross],
                      'scripts_called': [self.script],
                      'output_files': [self.results_filename],
                      'run_date': date_str
                      }
        with open(self.results_filename + '.json', 'w') as outfile:
            json.dump(self.myrec, outfile, indent=2, sort_keys=False)
            outfile.close()

    def cleanup(self):
        StatisticsFromCounts.remove_file(self, self.unsorted_filename)

    def run_all(self):
        rlogger = StatisticsFromCounts.logger
        rlogger.info("Script parameter= %s", self.script)
        rlogger.info("Cross parameter= %s", self.cross)
        # self.DATASET = self.cross
        self.BASENAME = os.path.basename(self.cross)
        self.script_run()
        self.load_raw_counts()
        self.merge_raw_stats()
        self.cleanup()
        self.run_recorder()

def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("script",
                        help="The R script [path/]<script>",
                        type=str)
    parser.add_argument("cross",
                        help="Data file prefix e.g. Col_x_Tsu",
                        type=str)
    parser.add_argument("--debug",
                        help="Adds messages to log and STDOUT.",
                        action="store_true")
    global args
    args = parser.parse_args()

# In case program is launched like this:
# $ python statistics_from_counts.py
if __name__ == "__main__":
    arg_parser()
    instance = StatisticsFromCounts(
                    args.script,
                    args.cross,
                    args.debug)
    instance.run_all()
