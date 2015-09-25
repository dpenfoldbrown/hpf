#!/usr/bin/env python

## Fragmentor class and associated functionality.
## Auth. Duncan Penfold-Brown, 1/18/2010.


## Imports
import os
import fileinput
import ConfigParser
from sys import stderr
from hpf.utilities import paths, copy, system


## The Fragmentor class, representing a sequence and fragmentation functionality and support for
#   that sequence.
# Contains a primary 'run' method responsible for running the entire fragment picking process. This
#   method calls support methods to initialize frag properties, set up working directory and file
#   structures, run the fragment picking toolchain, filter and package final results, and clean up
#   the system after completion (or failure).
class Fragmentor():

    def __init__(self, fasta_file="default.fasta", code="ZZ", results_dir="/scratch/fragmentor/results", config_file=None):
    # Initializes all variables and properties, except for those found in the _set_results
    # function (results instance variables) and the _parse_config function (tool and DB variables).
    # If config_file is not specified, attempts to use default config file 'fragmentor.config' in CWD.
        
        self.fasta_file = fasta_file
        self.code = code
        
        # Set name to fasta file name, without extension, via basename and string partitioning.
        self.name = os.path.basename(fasta_file).partition(".fasta")[0]
                
        # Set results directory and subdirectory: <results_dir>/<code>/<name>
        self.results_dir = os.path.join(os.path.expanduser(results_dir), self.code, self.name)
        
        # Set working directory for scratch files: <results_dir>/jobs/<name>
        self.working_dir = os.path.join(os.path.expanduser(results_dir), "jobs", self.name)
        
        # Set fragment picking results instance variables - blank, populated via self.run()
        self.frag3_lib = None
        self.frag9_lib = None
        self.frag_psipred  = None
        self.frag_psipred2 = None
        self.frag_fasta    = None
        self.results_files = None
        
        # If config file not specified, set to 'fragmentor.config' in the same dir as the module.
        if not config_file:
            self.config_file = os.path.dirname(__file__)+"/fragmentor.config"
        else:
            self.config_file = config_file
        
        # Parse config. file to get tool and database locations.
        self._parse_config(self.config_file)
        
    
    def run(self, nohoms=False, cleanup=True):
        
        try:
            # Populate environment with tool & DB locations (req'ed by make_fragments.local.pl).
            print "Fragmentor - run: set environment"
            self._set_env()
            
            # Create working and results directory structures, and move into working directory.
            print "Fragmentor - run:set directories"
            self._dir_setup()
            
            # Copy fasta file to working directory, set as working_fasta_file.
            print "Fragmentor - run: copy fasta files to working dir"
            copy(self.fasta_file, self.working_dir)
            working_fasta_file = os.path.join(self.working_dir, os.path.basename(self.fasta_file))
            
            # Run make_fragments script (treat as a black box).
            print "Fragmentor - run: run frag script"
            self._run_frag(nohoms, working_fasta_file)
            
            # Check for results files and set results instance variables.
            print "Fragmentor - run: set results variables and structure"
            self._set_results()
            
            # Filter 3' and 9' fragment library results (cut column size, bzip files).
            print "Fragmentor - run: filter fragment libs"
            self._filter_fraglibs()
            
            # Gather results (move all files to results directory).
            print "Fragmentor - run: gather results"
            self._gather_results()
        
        except:
            print "\nFragment picking failed. Exiting.\n"
            raise
        
        finally:    
            if cleanup:
                # Clean up work files and finish.
                self.cleanup()
        
        print "\nFragment picking completed successfully\n"
        
    
    def _set_env(self):
        # Put local variables rep'ing tools and DBs into the environment (via os.environ).
        env = {
                'BLAST_DIR': self.blast,
                'NNMAKE_SHORT_DIR': self.nnmake_short,
                'PSIPRED_DIR': self.psipred,
                'JUFO_DIR': self.jufo,
                'SAM_DIR': self.sam,
                'SAM_2ND_DIR': self.sam_2nd,
                'NR_DIR': self.nr,
                'NNMAKEDB_DIR': self.nnmakedb,
              }
        for key in env:
            # If paths exist, put in Environment (otherwise, fail).
            paths.existsOrFail(env[key])
            os.environ[key] = env[key]
    
    
    def _dir_setup(self): 
        # Create results and working directories.
        paths.ensure(self.results_dir)
        paths.ensure(self.working_dir)
        
        #try:
        #    os.makedirs(self.results_dir)
        #    os.makedirs(self.working_dir)
        #except OSError as e:
        #    print e.errno, e.strerror
        #    raise

        # Change directory to working directory.
        os.chdir(self.working_dir)
    
    
    def _run_frag(self, nohoms, fasta):
    # Takes a boolean nohoms (to pass to script) and the fasta input (filename).
    
        # Check to see if script given in config file exists.
        if not os.path.isfile(self.frag_script):
            raise Exception("Fragment picking {0} script specified in config file ({1}) is not" \
                            "accessible.".format(self.frag_script, cofig_file))
        
        # Run the given script on fasta file. Current options: ?-nohoms -nosam -verbose.
        if nohoms:
            cmd = "{0} -nohoms -nosam -verbose {1}".format(self.frag_script, fasta)
        else:
            cmd = "{0} -nosam -verbose {1}".format(self.frag_script, fasta)
        
        print "Fragmentor: executing frag command {0}".format(cmd)
        system(cmd)
    
    
    def _set_results(self):
    # All static functionality. If results files change (new fragment picking software), must 
    # manually recode (sorry).
    # Sets instance variables representing fragment picking results: frag3_lib, frag9_lib,
    # frag_psipred, frag_psipred2, frag_fasta, and results_files
    
        # Results from fragment picking run should be the following five files (code xx, number 001):
        # aaxx00103_05.075_v1_3, aaxx00109_05.075_v1_3 - 3' and 9' fragment libraries.
        # xx001.psipred, xx001.psipred_ss2 - psipred prediction files.
        # xx001.fasta - original fasta sequence file.
        
        self.frag3_lib = os.path.join(self.working_dir, "aa"+self.name+"03_05.075_v1_3")
        self.frag9_lib = os.path.join(self.working_dir, "aa"+self.name+"09_05.075_v1_3")
        self.frag_psipred  = os.path.join(self.working_dir, self.name+".psipred")
        self.frag_psipred2 = os.path.join(self.working_dir, self.name+".psipred_ss2")
        self.frag_fasta    = os.path.join(self.working_dir, self.name+".fasta")
        
        self.results_files = [self.frag3_lib, self.frag9_lib, self.frag_psipred, self.frag_psipred2, self.frag_fasta]
        
        for res in self.results_files:
            if not os.path.isfile(res):
                raise Exception("Results file {0} does not exist - problem with fragment picking run.".format(res))

    
    def _filter_fraglibs(self):
    # Runs on fragment picking result 3' and 9' fragment libraries (must exist from Frag run).
    # Filters the file inplace: writes " " + first 47 chars of the line. (IBM formatting).
    # Also uses bzip2 binary to compress fragment library files (IBM formatting). Can do with bz2 library.
        
        if self.frag3_lib == None or self.frag9_lib == None:
            raise Exception("Fragment library results files do not yet exist. Run Fragmentor.run(...) first.")
        
        frag_libs = [self.frag3_lib, self.frag9_lib]
        for frag_lib in frag_libs:
            # Filter frag library columns.
            for line in fileinput.input(frag_lib, inplace=1):
                print " "+line[:47].strip()
            fileinput.close()
            
            # Bzip frag library file and reset instance variable to reflect new name.
            cmd = "bzip2 {0}".format(frag_lib)
            system(cmd)
         
        # Update instance variable values to reflect bzip on fragment files.
        self.results_files.remove(self.frag3_lib)
        self.results_files.remove(self.frag9_lib)
        self.frag3_lib = self.frag3_lib+".bz2"
        self.frag9_lib = self.frag9_lib+".bz2"
        self.results_files.append(self.frag3_lib)
        self.results_files.append(self.frag9_lib)
    
    
    def _gather_results(self):
    # Move all results files to results directory and update instance variable values to reflect 
    # new location.
        
        for res in self.results_files:
            copy(res, self.results_dir)
        
        self.frag3_lib = os.path.join(self.results_dir, os.path.basename(self.frag3_lib))
        self.frag9_lib = os.path.join(self.results_dir, os.path.basename(self.frag9_lib))
        self.frag_psipred  = os.path.join(self.results_dir, os.path.basename(self.frag_psipred))
        self.frag_psipred2 = os.path.join(self.results_dir, os.path.basename(self.frag_psipred2))
        self.frag_fasta    = os.path.join(self.results_dir, os.path.basename(self.frag_fasta))
        
        self.results_files = [self.frag3_lib, self.frag9_lib, self.frag_psipred, self.frag_psipred2, self.frag_fasta]
        
    
    def cleanup(self):
        # Remove working directory and all its contents.
        paths.removerf(self.working_dir)


    def _parse_config(self, config_file):
        config = ConfigParser.ConfigParser()
        
        # Try to read given config file.
        try:
            if not os.path.isfile(config_file):
                raise IOError("File does not exist.")
            config.read(config_file)
            
        except IOError:
            print >> stderr, "Error in reading config file: {0}. Check that file exists and is readable.".format(config_file)
            raise
        except ConfigParser.ParsingError:
            print >> stderr, "Error parsing config file: {0}. Check for correct formatting.".format(config_file)
            raise
        except:
            print >> stderr, "Error with config file: {0}.".format(config_file)
            raise
            
        # Parse config parameters into local variables.
        try:
            self.frag_script  = os.path.expanduser(config.get("script", "frag_script"))
            self.blast        = os.path.expanduser(config.get("tools", "blast"))
            self.nnmake_short = os.path.expanduser(config.get("tools", "nnmake_short"))
            self.psipred      = os.path.expanduser(config.get("tools", "psipred"))
            self.jufo         = os.path.expanduser(config.get("tools", "jufo"))
            self.sam          = os.path.expanduser(config.get("tools", "sam"))
            self.sam_2nd      = os.path.expanduser(config.get("tools", "sam_2nd"))
            self.nr       = os.path.expanduser(config.get("databases", "nr"))
            self.nnmakedb = os.path.expanduser(config.get("databases", "nnmakedb"))
            
        except ConfigParser.InterpolationError:
            print >> stderr, "Error interpolating substituted values in config file: {0}. Check interpolated values.".format(config_file)
            raise
        except ConfigParser.NoSectionError, ConfigParser.NoOptionError:
            print >> stderr, "Section or Option does not exist. Check config file: {0}".format(config_file)
            raise
        except:
            print >> stderr, "Error getting values from config file. Check config file {0}".format(config_file)
            raise

