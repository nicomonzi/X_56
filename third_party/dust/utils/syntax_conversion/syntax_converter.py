#./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\.
#.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////..
#..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\......
#...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\......
#....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\......
#.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\......
#......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\......
#.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
#........\///////////........\////////......\/////////..........\///.......
##=========================================================================
##
## Copyright (C) 2018-2022 Politecnico di Milano,
##                           with support from A^3 from Airbus
##                    and  Davide   Montagnani,
##                         Matteo   Tugnoli,
##                         Federico Fonte
##
## This file is part of DUST, an aerodynamic solver for complex
## configurations.
##
## Permission is hereby granted, free of charge, to any person
## obtaining a copy of this software and associated documentation
## files (the "Software"), to deal in the Software without
## restriction, including without limitation the rights to use,
## copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the
## Software is furnished to do so, subject to the following
## conditions:
##
## The above copyright notice and this permission notice shall be
## included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
## EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
## OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
## HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
## WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
## OTHER DEALINGS IN THE SOFTWARE.
##
## Author: Andrea Colli
##=========================================================================

# Script for DUST syntax conversion
# Run with -h to display help text

import os
import sys
import argparse
from pathlib import Path
import shutil
import csv
import re

info_text = 'This script converts the syntax used in DUST files. It can be used on input (.in) or source (.f90) files.\
            Starting from the specified folder, it searches recursively all subfolders for the files with the appropriate extension\
            and it replaces all the strings in each file as described in the input dictionary file. The dictionary file is\
            a CSV file containing for each row a couple of strings in the format \'new_string,old_string\'; the script\
            replaces every occurrence of \'old_string\' with \'new_string\'.\n WARNING! USE WITH CAUTION: the script operates\
            on ALL files and folders with the specified extension or name, regardless whether they are actually DUST-related.\
            It is the user\'s responsibility to check for unwanted behaviour; for this, a dry run option is provided.'

# Parse command line arguments
parser = argparse.ArgumentParser(description='converts DUST input files syntax according to the key pairs in dictionary file',\
                                epilog=info_text)
parser.add_argument('--delete-backup', action="store_true", help='delete backup folders from previous runs')
parser.add_argument('--no-backup', action="store_true", help='do not backup files')
parser.add_argument('-i', '--invert', action="store_true", help='convert back from new syntax to old syntax')
parser.add_argument('-s', '--source', action="store_true", help='apply conversion to source files; overrides -p and -b')
parser.add_argument('-n', '--dry-run', action="store_true", help='perform dry run test (without modifying files)')
parser.add_argument('-p', '--pattern', help='pattern used to find files with glob; use with caution', default='*.in')
parser.add_argument('-b', '--backup-dirname', help='name for the backup directory; use with caution', default='old_in_files')
parser.add_argument('-d', '--dictionary', help='file containing key pairs')
parser.add_argument('folder', help='root folder from which to start the recursive search')
args = parser.parse_args()

# Check if script is run as root, abort for safety
if os.geteuid() == 0:
    print('ERROR: running script with root privileges, this is not allowed for safety reasons.')
    print('Exiting...')
    exit()     

# Dry run
if args.dry_run:
    print('Executing dry run, no file will be modified.\n')

# Root folder
base_folder = args.folder
# Pattern for finding files
pattern = args.pattern

# Name of the backup directory to be created
backup_name = args.backup_dirname
if args.source: # Override for --source switch
    pattern = '*.f90'
    backup_name = 'old_src_files'

print(backup_name)
if args.delete_backup: # Delete backup folders
    # Ask for confirmation
    print('WARNING: all backups in this folder and all subfolders will be deleted.')
    flag = input('Continue? [y/N] ')
    if not(str(flag).lower() == 'y' or str(flag).lower() == 'yes'):
        print('Exiting...')
        exit() # Stop execution if not confirmed
    else:
        for bak_dir in Path(base_folder).rglob(backup_name):
            if bak_dir.is_dir():
                print('Deleting '+str(bak_dir))
                if not args.dry_run:
                    shutil.rmtree(bak_dir)
    print('\nFinished. Exiting...')
    exit()   

elif args.dictionary is None: # Check if dictionary file is specified
    print('ERROR: dictionary file not specified. Exiting...')
    exit()
    
else: # Perform substitution
    # CSV dictionary file with the keys to be replaced, format:
    # old_key,new_key    
    dict_name = args.dictionary
    print('Reading dictionary file...\n')
    # Import key pairs from file
    with open(dict_name, 'r') as f:
        keys = [tuple(line) for line in csv.reader(f)]

    if args.no_backup: # Ask for confirmation
        print('WARNING: no backup of the files will be created.')
        flag = input('Continue? [y/N] ')
        if not(str(flag).lower() == 'y' or str(flag).lower() == 'yes'):
            print('Exiting...')
            exit() # Stop execution if not confirmed
    
    # Recursively find .in files in base folder and subfolders and perform substitution
    print('Performing substitution...\n')
    for file in Path(base_folder).rglob(pattern):    
            
        if file.resolve() is Path(sys.argv[0]).resolve(): # Avoid having the script modify itself
            print('ERROR: the script is attempting to modify itself, this is not allowed for safety reasons.')
            print('Exiting...')
            exit()            
        if '.bak' in str(file):
            print('ERROR: the script is attempting to modify a probable backup file.')
            print('Exiting...')
            exit()
            
        print('Processing file '+str(file))
        pwd = file.parent        
        
        if not args.no_backup: # Create backup
            # Create backup folder
            backup_folder = pwd.resolve().joinpath(backup_name)
            print('Saved old files in '+str(backup_folder)+'\n')
            if not args.dry_run:
                backup_folder.mkdir(exist_ok=True)
            # Copy file to backup folder and rename it
            backup_file = backup_folder.joinpath(str(file.name)+'.bak')
            if not args.dry_run:
                shutil.copy(file,backup_file)     

        # Read input file and replace keys, store into string list
        content = []
        with open(file,'r') as f:
            for line in f.readlines():
                for old, new in keys:
                    if args.source:
                        regexp_old = re.compile(r"(?i)'\b"+old+r"\b'")
                        regexp_new = re.compile(r"(?i)'\b"+new+r"\b'")
                        old = '\''+old+'\''
                        new = '\''+new+'\''
                    else:
                        regexp_old = re.compile(r"(?i)\b"+old+r"\b")
                        regexp_new = re.compile(r"(?i)\b"+new+r"\b")
                    # only change line for .in file or in CreateOption for source file
                    modify_source_line = args.source and ('%Create' in line or '=get' in line.replace(' ','') or\
                    'countoption' in line or 'check_opt_consistency' in line)
                    if not args.source or (args.source and modify_source_line):
                        if args.invert: # inverse behaviour, replace new key with old
                            #line = line.replace(new, old)
                            line = re.sub(regexp_new, old, line)
                        else: # normal behaviour, replace old key with new
                            line = re.sub(regexp_old, new, line)
                            #line = line.replace(old, new)
                content.append(line)
                            
        # Overwrite input file
        if not args.dry_run:
            with open(file, 'w') as f:
                for line in content:
                    f.writelines(line)
    
    print('Finished. Exiting...')
    
