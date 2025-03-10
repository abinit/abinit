#!/usr/bin/env python

#
# Copyright (C) 2010-2025 ABINIT Group (Jean-Michel Beuken)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# 'false' errors elimination in the 'linkchecker_ext.log' generated by
# a script located on ref slave ( abiref:~buildbot/bin/LinkChecker.sh )
#

from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os
import re
from lxml import etree
import requests
import argparse

# ---------------------------------------------------------------------------- #

#all_tags = [ 'url','name','parent','realurl','extern','dlsize','checktime','level','infos','valid' ]
#printable_tags = [ 'url', 'name', 'parent', 'realurl', 'valid' ]

version_info = (1, 0, 0)
version = '.'.join(str(c) for c in version_info)

debug = False
server = "http://localhost:8000"

#
# Functions
#

def rm_server(keyword):
  if keyword.startswith(server):
     return keyword[len(server):]

def Checking_on_url_to_skip(e, u, v):
  global url_to_skip
  for ui in url_to_skip :
     if e == "1" and ui == u and  ( v == "syntax OK" or v == "filtered" ) :
         #print("url_to_skip")
         return True
  return False


def Checking_on_url_string_to_skip(e, u):
  global url_string_to_skip
  for s in url_string_to_skip :
     if e == "1" and u.find( s ) >= 0 :
         #print(s,"url_string_to_skip")
         return True
  return False


def Checking_on_no_error_list(url, info, valid):
  global no_error_list
  #print("Enter no_error_list... %s,%s,%s" % (url.text,info,valid))
  for no_error in no_error_list:

    norc = True  # True : we consider that this is not an error

    url_rc = None
    try:
       url_rc = no_error['url'].search(url.text)
    except:
       print("-- no URL NAME --")
    norc = norc and ( url_rc != None )

    info_rc = None
    if info is not None:
      try:
       info_rc = no_error['info'].search(info.text)
       norc = norc and ( info_rc != None )
      except:
       pass
    else:
       info_rc = True

    valid_rc = None
    valid_rc = no_error['valid'].search(valid)
    norc = norc and ( valid_rc != None )

    if url_rc != None and info_rc != None and valid_rc != None:
        return True # in the exception list -> next xml entry
  return False # may be a error


def Checking_on_false_error_list(url, valid, parent, cnx):
  global false_error_list

  try:
     URL=url.text
  except:
     return False

  for false_error in false_error_list:

    url_rc = None
    try:
       url_rc = false_error['url'].search(URL)
    except:
       url_rc = True

    valid_rc = None
    try:
       valid_rc = false_error['valid'].search(valid)
    except:
       valid_rc = True

    parent_rc = None
    try:
       parent_rc = false_error['parent'].search(parent)
    except:
       parent_rc = True

    if cnx == 0:
       if url_rc != None and valid_rc != None and parent_rc != None:
          return True # is a false error
    else:
       if url_rc != None and valid_rc != None and parent_rc != None and cnx == 404:
          #print("found a false error...")
          return True # is a false error

  return False # may be a error


def Checking_on_warning_error_list(valid) :
  global warning_list

  for warning_error in warning_list:
    
    v_rc = None
    v_rc = warning_error['valid'].search(valid) 
    if v_rc != None :
       return True

  return False

# ---------------------------------------------------------------------------- #

# Types of error :
#    Error: 401 Unauthorized
#    Error: 403 Forbidden
#    Error: 404 Not Found
#    Error: 504 Gateway Time-out
#    Error: 502 Bad Gateway
#    Error: ReadTimeout:
#    ConnectionError: ('Connection aborted.

warning_list = [
    { 'valid': re.compile('^ReadTimeout')
    },
]

no_error_list = [
    { 'url'  : re.compile('(doi|aps|stacks.iop).org'),
      'info' : re.compile('^Redirected'),
      'valid': re.compile('^403 Forbidden')
    },
    { 'url'  : re.compile('jstor.org'),
      'valid': re.compile('^403 Unauthorized')
    },
]

false_error_list = [
    { 'url'  : re.compile('(dx.doi.orgg|en.wwikipedia.org)'),
      'valid': re.compile('^ConnectionError'),
    },
    { 'url'  : re.compile('10.1102/physrevb.27.4760'),
      'cnx'  : 404
    },
    { 'url'  : re.compile('10.1103/physrevb.87.085323'),
      'cnx'  : 404
    },
    { 'url'   : re.compile('abiconfigg'),
      'parent': re.compile('testlink/'),
      'cnx'   : 404
    },
    { 'url'   : re.compile('FAKE_URL'),
      'parent': re.compile('testlink/')
    },
]

url_to_skip = [
    "https://github.com/abinit/abiconfig",
    "https://github.com/abinit/abiflows",
    "https://github.com/abinit/abiconda",
    "https://github.com/abinit/abiconfig",
    "https://github.com/abinit/abitutorials",
    "https://github.com/abinit/abipy",
    "https://github.com/abinit/abiout",
    "https://github.com/abinit/abinit/",
    "https://github.com/abinit/abinit",
    "https://github.com/abinit/",
    "https://github.com/abinit",
    "https://www.facebook.com/abinit.org",
    "https://fonts.gstatic.com"
]

url_string_to_skip = [
    "cdn.jsdelivr.net",
    "cdn.embedly.com",
    "cdn.plot.ly",
    "maxcdn.bootstrapcdn.com",
    "facebook.com/abinit",
    "github.com/abinit/abinit/edit",
    "markdown-here/wiki/Markdown",
    "nschloe/betterbib",
    "github.com/mitya57",
    "github.com/helderco",
    "github.com/abinit/abinit/commit/",
    "github.com/abinit/abipy_assets",
    "github.com/abinit/abinit/tree",
    "github.com/abinit/abinit/blob",

]

# ---------------------------------------------------------------------------- #

#
# Main program
#
def main(filename,home_dir=""):
  from os.path import join as pj

  # Check if we are in the top of the ABINIT source tree
  my_name = os.path.basename(__file__) + ".main"
  if ( not os.path.exists(pj(home_dir,"configure.ac")) or
       not os.path.exists(pj(home_dir, "src/98_main/abinit.F90")) ):
    print("%s: You must be in the top of an ABINIT source tree." % my_name)
    print("%s: Aborting now." % my_name)
    sys.exit(1)

  #
  tree = etree.parse(filename)
  
  rc=0 # true error counter 
  frc=0 # false error counter
  wrc=0 # warning error counter

  urls=set()

  for child in tree.xpath("/linkchecker/urldata"):
  
    url    = child.find('url')
    parent = child.find('parent')
    URL    = url.text
    info   = child.find('infos/info')
    extern = child.find('extern')
    valid  = child.find('valid').get("result")

    ### check for duplicate entry except for FAKE_URL

    try:
      if not ("FAKE_URL" in URL) :
          if URL in urls :
             continue
          else: 
             urls.add(URL)
      #else:
      #    if not ("index.html" in parent.text) :
      #       continue
    except:
      pass

    ### precleaning ###

    v = re.compile("^200")   # status "200" or "200 OK"
    if v.search(valid):
       continue

    if url.text[0:6] == "mailto" and valid == "Valid mail address syntax" :
       continue

    if Checking_on_url_to_skip( extern.text, url.text, valid ):
       continue

    if Checking_on_url_string_to_skip( extern.text, url.text ):
       continue

    ### fine cleaning ###

    if Checking_on_no_error_list(url, info, valid) :
       continue

    ### last chance to know if it's not a error ###
    ### access denied  then checks with 'curl'  ###

    Check_connection = False
    cnx_status = 0
    if valid == "syntax OK" :
        Check_connection = True
        if debug : print("check cnx : ",url.text)
        try: 
           request = requests.get(url.text, headers={"content-type":"text"}, timeout=(2,2) )
           cnx_status = request.status_code
        except (requests.Timeout, requests.ConnectionError, KeyError) as e:
           if debug : print('failed to connect to website ({})'.format(e))
           continue
        if cnx_status == 200 :  # OK
           continue
        if cnx_status == 403 :  # cnx ok but Forbidden for robot
           continue
    
    # check if the error is a "false" error
    if Checking_on_false_error_list(url=url, valid=valid, parent=parent.text, cnx=cnx_status) :
        frc += 1
        continue

    # check if it's a warning
    if Checking_on_warning_error_list(valid=valid) :
        wrc += 1

    # found a true error... : reporting on bb
    rc += 1
    name=child.find('name')
    realurl=child.find('realurl')
    try:
       print("{0:12} {1}".format('URL',url.text))
    except:
       print("{0:12} {1}".format('URL',' ** NO URL **'))
    try:
       print("{0:12} {1}".format('Name',name.text))
    except:
       print("{0:12} {1}".format('Name','** NO NAME **'))
    print("{0:12} {1}, line {2}".format('Parent URL',rm_server(parent.text),parent.get('line')))
    try:
       print("{0:12} {1}".format('Infos',info.text))
    except:
       pass 
    print("{0:12} {1}".format('Real URL',realurl.text))
    print("{0:12} {1}".format('Result',valid))
    if Check_connection : 
          print("{0:12} {1}".format('Status CNX',request.status_code))
    print('---------------------------')

  print('false errors found : ',frc,' ( must be : 7 )')

  if rc != 0:
     if rc - wrc != 0:
        return 2 # FAILED
     else:
        print('warning errors found : ',wrc, ' [probably, these errors are transient...]')
        return 1 # WARNING

  return 0 # SUCCESS

# ---------------------------------------------------------------------------- #

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Remove false errors')
  parser.add_argument('--verbose', '-v', action='count',
                      help='increase verbosity. Specify multiple times')
  parser.add_argument('--version', action='version',
                      version='%(prog)s {}'.format(version),
                      help='show the version number and exit')
  parser.add_argument('filename', help='input file (xml format)'),
  parser.add_argument('home_dir', nargs='?', default=os.getcwd())

  args = parser.parse_args()

  filename = args.filename
  home_dir = args.home_dir
  #print(filename, home_dir)

  exit_status = main(filename,home_dir)
  sys.exit(exit_status)
