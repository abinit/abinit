#!/usr/bin/python3

import os
import re
import sys

if ( sys.version_info[0] == 3 ):
    from configparser import ConfigParser
else:
    from ConfigParser import ConfigParser

                    ########################################

#
# Case-sensitive configuration parser
#

class AbinitConfigParser(ConfigParser):

    def __init__(self, config_text):

        # Make parent happy first
        ConfigParser.__init__(self)

        # Read configuration from file or text string
        if ( (not re.search("\n", config_text, flags=re.DOTALL)) and \
             os.access(config_text, os.R_OK) ):
            self.read(config_text)
        else:
            if ( sys.version_info[0] == 3 ):
                self.read_string(config_text)
            else:
                raise TypeError("use Python 3 to read configurations from strings")


    def optionxform(self, option):
    
        return str(option)


                    ########################################


#
# Configuration extractor
#

class AbinitConfigExtractor(object):

    def __init__(self, prefix, config_file):

        self.parser = AbinitConfigParser(config_file)


                    ########################################

#
# Shell-script generator
#

class AbinitConfigShellMaker(object):

    def __init__(self, vendor_cases):

        if ( not isinstance(vendor_cases, dict) ):
            raise TypeError("vendor_cases must be a dictionary")

