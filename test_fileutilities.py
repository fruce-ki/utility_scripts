#!/usr/bin/python

"""Unittests for fileutilities.py

Author: Kimon Froussios
Compatibility tested: python 2.7.3
Last reviewed: 12/04/2016

The tests generate many temporary files with names such as 'testfile' and 
variations and composites thereof in the current directory. Please ensure you 
have no such named files that you wish to keep, as they will be overwritten and 
subsequently deleted.

The tests are by no means exhaustively complete. Priority was given to
verifying that things behave as expected under normal operation parameters and
even then the tests are not exhaustive. But all of them combined, they should 
directly or indirectly cover everything, including the embedded main() of 
fileutilities.py. Exceptional circumstances and abnormal parameters are not
tested, although other than bad input data format there shouldn't be many
things that can normally go wrong.
"""

import unittest, os, sys, zlib, subprocess, re
import fileutilities as fu
import pandas as pd


# Frequently used value:
_invalidPath = "really/long/and/highly/improbable/testfilename"
scriptdir = os.path.split(sys.argv[0])[0]
script = os.path.join(scriptdir, "fileutilities.py")


  
    
class Test_expandPaths(unittest.TestCase):
    def test_expand(self):
        """Test path expansion for mutliple files."""
        # Get some reference values
        validpath = os.path.abspath(os.path.expanduser(sys.argv[0]))
        invalidpath = os.path.abspath(os.path.expanduser(_invalidPath))
        homepath = os.path.abspath(os.path.expanduser("~/"))
        # Expand a valid path, an invalid path and one that needs expansion:
        self.assertEqual(fu.expand_fpaths([sys.argv[0], _invalidPath, "~"]), 
                        [validpath, invalidpath, homepath])
        
      
class Test_FilesListBasics(unittest.TestCase):
    def setUp(self):
        # Prepare some files. They don't need to actually exist.
        self.files = ["testfile1", "testfile2", "testfile3"]
        self.aliases = ["prefx_suffx", "alias", "testfile3"]
        # FileUtilities has a function for this, but since this is a test module, do it the hard way.
        self.expanded = [os.path.abspath(os.path.expanduser(self.files[0])),
                         os.path.abspath(os.path.expanduser(self.files[1])),
                         os.path.abspath(os.path.expanduser(self.files[2]))]    
    # Prepare a list with aliases
        self.listfile2 = "listtestfile2"
        with open(self.listfile2, 'w') as f:
            f.write(self.files[0] + "\tprefx\tsuffx\n" +
                    self.files[1] + "\talias\n" +
                    self.files[2] + "\n")
          
    def tearDown(self):
        os.remove(self.listfile2)
          
    def test_instantiationDefault(self):
        """Test the creation of empty instance."""
        flist = fu.FilesList()
        self.assertEqual(flist[0:3], [])
        self.assertEqual(flist.aliases[0:3], [])
          
    def test_instantiationFromList(self):
        """Test the conversion of a built-in list to a FilesList."""
        flist = fu.FilesList(self.files)
        self.assertEqual(flist[0:3], self.expanded)
        self.assertEqual(flist.aliases[0:3], self.files)
      
    def test_instantiationFromFilesAndAliases(self):
        """Test creating a FilesList from a list of files and a list of aliases"""
        flist = fu.FilesList(self.files, self.aliases)
        self.assertEqual(flist[0:3], self.expanded)
        self.assertEqual(flist.aliases[0:3], self.aliases)
      
    def test_instantiation_fromTuples(self):
        flist = fu.FilesList(fromtuples=zip(self.files,self.aliases))
        self.assertEqual(flist[0:3], self.expanded)
        self.assertEqual(flist.aliases[0:3], self.aliases)
        # Test tuples priority by swapping the values for flist and aliases (both should be ignored in favour of the tuples data).
        flist = fu.FilesList(self.aliases, self.files, fromtuples=zip(self.files,self.aliases))
        self.assertEqual(flist[0:3], self.expanded)
        self.assertEqual(flist.aliases[0:3], self.aliases)
              
    def test_appendFile(self):
        """Test appending with one argument, like a list."""
        flist = fu.FilesList(self.files)
        flist.append("somefile")
        self.assertEqual(flist[-1], os.path.abspath("somefile"))
        self.assertEqual(flist.aliases[-1], "somefile")
              
    def test_appendFileWithAlias(self):
        """Test appending path and alias to the FilesList."""
        flist = fu.FilesList(self.files)
        flist.append("somefile", "rubbish")
        self.assertEqual(flist[-1], os.path.abspath("somefile"))
        self.assertEqual(flist.aliases[-1], "rubbish")
          
    def test_get(self):
        """Test representation of (path,alias) pairs as tuples."""
        flist = fu.FilesList(self.files)
        self.assertEqual(flist.get(2), (self.expanded[2], self.files[2]))
              
    def test_enumeration(self):
        """Test the custom enumerator."""
        flist = fu.FilesList(self.files)
        for f, (myfile, myalias) in flist.enum():
            self.assertEqual(myfile, self.expanded[f])
            self.assertEqual(myalias, self.files[f])
          
    def test_stringRepresentation(self):
        """Test the representation of an FilesList as a string."""
        flist = fu.FilesList().populate_from_files([self.listfile2])
        s = "0\t" + self.expanded[0] + "\t" + self.aliases[0] + "\n1\t" + self.expanded[1] + "\t" + self.aliases[1] + "\n2\t" + self.expanded[2] + "\t" + self.aliases[2] + "\n"
        self.assertEqual( str(flist), s)
           
      
class Test_FilesListFileAccess(unittest.TestCase):
    def setUp(self):
        # Prepare some files. They don't need to actually exist.
        self.files = ["testfile1", "testfile2", "testfile3"]
        self.aliases = ["prefx_suffx", "alias", "testfile3"]
        self.aliases2 = ["prefx_suffx_2", "alias_2", "testfile3_2"]
        # FileUtilities has a function for this, but since this is a test module, do it the hard way.
        self.expanded = [os.path.abspath(os.path.expanduser(self.files[0])),
                         os.path.abspath(os.path.expanduser(self.files[1])),
                         os.path.abspath(os.path.expanduser(self.files[2]))]
        # Prepare a plain list of files
        self.listfile = "listtestfile"
        with open(self.listfile, 'w') as f:
            f.write("\n".join(self.files))
        # Prepare a list with aliases
        self.listfile2 = "listtestfile2"
        with open(self.listfile2, 'w') as f:
            f.write(self.files[0] + "\tprefx\tsuffx\n" +
                    self.files[1] + "\talias\n" +
                    self.files[2] + "\n")
        self.listfile3 = "listtestfile3"
        # And one with a different field separator.
        with open(self.listfile3, 'w') as f:
            f.write(self.files[0] + ",prefx,suffx\n" +
                    self.files[1] + ",alias\n" +
                    self.files[2] + "\n")
        # Placeholder for output file
        self.outfile = "outlisttestfile"
        with open(self.outfile, 'w'):
            pass
          
    def tearDown(self):
        os.remove(self.listfile)
        os.remove(self.listfile2)
        os.remove(self.listfile3)
        os.remove(self.outfile)
              
    def test_inputFromFileSimple(self):
        """Input file is a simple list of files."""
        flist = fu.FilesList().populate_from_files([self.listfile])
        self.assertEqual(flist[0:3], self.expanded)
        self.assertEqual(flist.aliases[0:3], self.files)
          
    def test_inputFromFileExtra(self):
        """Input file contains aliases in mixed format."""
        flist = fu.FilesList().populate_from_files([self.listfile2])
        self.assertEqual(flist[0:3], self.expanded)
        self.assertEqual(flist.aliases[0:3], self.aliases)
        # Append more files, and try different column separator.
        flist.populate_from_files([self.listfile3], colSep=",")
        self.assertEqual(flist[0:3], self.expanded)
        self.assertEqual(flist[3:6], self.expanded)
        self.assertEqual(flist.aliases[0:3], self.aliases2)
        self.assertEqual(flist.aliases[3:6], self.aliases)        
          
    def test_outputToFileAndBack(self):
        """Test outputting a FilesList to text file in format compatible for input."""
        flist  = fu.FilesList(self.files)
        flist.to_file(self.outfile)
        flist2 = fu.FilesList().populate_from_files([self.outfile])
        self.assertEqual(flist, flist2)
              
      
class Test_FilesListDirectoryAccess(unittest.TestCase):
    def setUp(self):
        self.files = ["testfile1.txt", "testingfile2", "testfile3a.pdf"]
        # Input directories
        self.d1 = "dirtestdirectory1"
        self.d2 = "dirtestdirectory2"
        self.d = os.getcwd()
        os.mkdir(self.d1)
        os.mkdir(self.d2)
        os.chdir(self.d1)
        with open(self.files[0], 'w'):
            pass
        with open(self.files[1], 'w'):
            pass
        os.chdir(self.d)
        os.path.exists(self.d2)
        os.chdir(self.d2)
        with open(self.files[2], 'w'):
            pass
        os.chdir(self.d)
          
    def tearDown(self):
        os.chdir(self.d1)
        os.remove(self.files[0])
        os.remove(self.files[1])
        os.chdir(self.d)
        os.rmdir(self.d1)
        os.chdir(self.d2)
        os.remove(self.files[2])
        os.chdir(self.d)
        os.rmdir(self.d2)
                  
    def test_filterFromDirectoriesDefaults(self):
        flist = fu.FilesList().populate_from_directories([self.d2, self.d1])
        self.assertEqual(len(flist), len(self.files))
        self.assertEqual(len(flist.aliases), len(self.files))
        for f in self.files:
            self.assertTrue(os.path.splitext(f)[0] in flist.aliases)
              
    def test_filterFromDirectoriesPattern(self):
        flist = fu.FilesList().populate_from_directories([self.d2, self.d1], patterns=["txt","ing"])
        self.assertTrue(os.path.splitext(self.files[0])[0] in flist.aliases)
        self.assertTrue(os.path.splitext(self.files[1])[0] in flist.aliases)
        flist = fu.FilesList().populate_from_directories([self.d2, self.d1], ["\d\."])
        self.assertTrue(os.path.splitext(self.files[0])[0] in flist.aliases)
      
      
class Test_TailOfFiles(unittest.TestCase):
    def setUp(self):
        self.testfile = "tailtestfile"
        with open(self.testfile, 'w') as fout:
            fout.write("line 1\nline 2\nline 3\n")
             
    def tearDown(self):
        os.remove(self.testfile)
                 
    def test_tailSmallerThanFile(self):
        """Test getting fewer lines than there are (including none)."""
        self.assertEqual(fu.tails([self.testfile], 0)[0], [])
        self.assertEqual(fu.tails([self.testfile], 1)[0], ["line 3\n"])
                     
    def test_tailEqualToFile(self):
        """Test getting all the lines."""
        self.assertEqual(fu.tails([self.testfile], 3)[0], ["line 1\n", "line 2\n", "line 3\n"])
                     
    def test_tailMoreThanFile(self):
        """Test getting more lines than there are."""
        self.assertEqual(fu.tails([self.testfile], 5)[0], ["line 1\n", "line 2\n", "line 3\n"])
             
    def test_tailTwoFiles(self):
        """Test tailing two (identical) files, to check output structure."""
        res = fu.tails([self.testfile, self.testfile], 2)
        self.assertEqual(res[0], res[1])     
             
    def test_tailInvalidPathRaisesException(self):
        """Test that exception is raised of a path is invalid."""
        self.assertRaises(IOError, fu.tails, [_invalidPath])        
                  
       
class Test_HeadlOfFiles(unittest.TestCase):
    def setUp(self):
        self.testfile = "headtestfile"
        with open(self.testfile, 'w') as fout:
            fout.write("line 1\nline 2\nline 3\n")
             
    def tearDown(self):
        os.remove(self.testfile)
                 
    def test_headSmallerThanFile(self):
        """Test getting fewer lines than there are (including none)."""
        self.assertEqual(fu.heads([self.testfile], 0)[0], [])
        self.assertEqual(fu.heads([self.testfile], 1)[0], ["line 1\n"])
                     
    def test_headEqualToFile(self):
        """Test getting all the lines."""
        self.assertEqual(fu.heads([self.testfile], 3)[0], ["line 1\n", "line 2\n", "line 3\n"])
                     
    def test_headMoreThanFile(self):
        """Test getting more lines than there are."""
        self.assertEqual(fu.heads([self.testfile], 5)[0], ["line 1\n", "line 2\n", "line 3\n"])
             
    def test_headTwoFiles(self):
        """Test tailing two (identical) files, to check output structure."""
        res = fu.heads([self.testfile, self.testfile], 2)
        self.assertEqual(res[0], res[1])     
             
    def test_headInvalidPathRaisesException(self):
        """Test that exception is raised of a path is invalid."""
        self.assertRaises(IOError, fu.heads, [_invalidPath])     
              
      
class Test_batchSymlink(unittest.TestCase):
    def setUp(self):
        self.fromdir = "sltestdirfrom"
        self.todir = "sltestdirto"    
        self.files = ["sltestfile1.a", "sltestfile2.b", "sltestfile3.a", "sltestfile4.b"]
        self.aliases = ["sltestlink1", "sltestlink2", "sltestlink3", "sltestlink4"]
        self.links = ["sltestlink1.a", "sltestlink2.b", "sltestlink3.a", "sltestlink4.b"]
        self.paths = [self.files[0], 
                      self.files[1], 
                      os.path.join(self.fromdir, self.files[2]),
                      os.path.join(self.fromdir, self.files[3])]
        self.d = os.getcwd()
        with open(self.files[0], 'w') as f:
            pass
        with open(self.files[1], 'w') as f:
            pass
        os.mkdir(self.fromdir)
        os.chdir(self.fromdir)
        with open(self.files[2], 'w') as f:
            pass
        with open(self.files[3], 'w') as f:
            pass
        os.chdir(self.d)
        os.mkdir(self.todir)
             
    def tearDown(self):
        # Current directory
        for i in [0, 1, 2, 3]:
            try:
                os.remove(self.files[i])
            except OSError:
                pass
        for i in [0, 1, 2, 3]:
            try:
                os.remove(self.aliases[i])
            except OSError:
                pass
        # Departure directory.
        os.chdir(self.fromdir)
        os.remove(self.files[2])
        os.remove(self.files[3])
        os.chdir(self.d)
        os.rmdir(self.fromdir)
        # Destination directory.
        os.chdir(self.todir)
        for i in [0, 1, 2, 3]:
            try:
                os.remove(self.files[i])
            except OSError:
                pass
        for i in range(0, 4):
            try:
                os.remove(self.links[i])
            except OSError:
                pass
        os.chdir(self.d)
        os.rmdir(self.todir)
        for i in range(0, 4):
            try:
                os.remove(self.links[i])
            except OSError:
                pass
          
    def test_symlinkFromFilesList(self):
        flist = fu.FilesList(files=self.paths, aliases=self.aliases)
        fu.slink(flist)
        self.assertTrue(os.path.exists(self.links[0]))
        self.assertTrue(os.path.exists(self.links[1]))
        self.assertTrue(os.path.exists(self.links[2]))
        self.assertTrue(os.path.exists(self.links[3]))
           
    def test_symlinkFromListAndAliases(self):
        fu.slink(self.paths, aliases=self.aliases)
        self.assertTrue(os.path.exists(self.links[0]))
        self.assertTrue(os.path.exists(self.links[1]))
        self.assertTrue(os.path.exists(self.links[2]))
        self.assertTrue(os.path.exists(self.links[3]))
           
    def test_symlinkWithoutAliases(self):
        # Can't use paths[0] and [1] as the default-named links in the default directory
        # would coincide with the actual files.
        fu.slink(self.paths[2:])
        self.assertTrue(os.path.exists(self.files[2]))
        self.assertTrue(os.path.exists(self.files[3]))
               
    def test_symlinkWithoutExtensions(self):
        # The following should already have extensions since they are named after the files.
        # Only using files[2] and [3] again since they are not in the target directory.
        fu.slink(self.paths[2:], autoext=False)
        self.assertTrue(os.path.exists(self.files[2]))
        self.assertTrue(os.path.exists(self.files[3]))
        # The following should not have extensions as they are named after the aliases.
        # I can use all files since the aliases won't overlap with the existing files.
        fu.slink(self.paths, aliases=self.aliases, autoext=False)
        self.assertTrue(os.path.exists(self.aliases[0]))
        self.assertTrue(os.path.exists(self.aliases[1]))
        self.assertTrue(os.path.exists(self.aliases[2]))
        self.assertTrue(os.path.exists(self.aliases[3]))
           
    def test_symlinkIntoDir(self):
        flist = fu.FilesList(files=self.paths, aliases=self.aliases)
        fu.slink(flist, dir=self.todir)
        self.assertTrue(os.path.exists(os.path.join(self.todir, self.links[0])))
        self.assertTrue(os.path.exists(os.path.join(self.todir, self.links[1])))
        self.assertTrue(os.path.exists(os.path.join(self.todir, self.links[2])))
        self.assertTrue(os.path.exists(os.path.join(self.todir, self.links[3])))
            
    def test_symlinkMainNoAlias(self):
        """Test command-line usage without aliases."""
        with open(os.devnull, 'w') as FNULL:
            command = ["python", script, 
                        "T", self.paths[0], self.paths[1], self.paths[2], self.paths[3],
                        "-C", "--link", self.todir ]
            subprocess.call(command, stdout=FNULL, shell=False)
        self.assertTrue(os.path.exists(os.path.join(self.todir, self.files[0])))
        self.assertTrue(os.path.exists(os.path.join(self.todir, self.files[1])))
        self.assertTrue(os.path.exists(os.path.join(self.todir, self.files[2])))
        self.assertTrue(os.path.exists(os.path.join(self.todir, self.files[3])))
                  
    def test_symlinkMainAliases(self):
        """Test command-line usage with aliases."""
        with open(os.devnull, 'w') as FNULL:
            command = ["python", script, 
                        "T", self.paths[0], self.paths[1], self.paths[2], self.paths[3], 
                        "-C", "--link", self.todir, self.aliases[0], self.aliases[1], self.aliases[2], self.aliases[3] ]
            subprocess.call(command, stdout=FNULL, shell=False)
        self.assertTrue(os.path.exists(os.path.join(self.todir, self.links[0])))
        self.assertTrue(os.path.exists(os.path.join(self.todir, self.links[1])))
        self.assertTrue(os.path.exists(os.path.join(self.todir, self.links[2])))
        self.assertTrue(os.path.exists(os.path.join(self.todir, self.links[3])))
                  
                   
class Test_LoopingArbitraryCommand(unittest.TestCase):
    def setUp(self):
        self.files = ["testdofile1", "testdofile2", "testdofile3.abc"]
        with open(self.files[0], 'w') as f:
            pass
        with open(self.files[1], 'w') as f:
            pass
        with open(self.files[2], 'w') as f:
            pass
           
    def tearDown(self):
        try:
            os.remove(self.files[0] + ".txt")
            os.remove(self.files[1] + ".txt")
            os.remove(self.files[2] + ".txt")
        except OSError:
            pass
        try:
            os.remove(self.files[0])
            os.remove(self.files[1])
            os.remove(self.files[2])
        except OSError:
            pass
        try:
            os.remove("testdofile3.txt")
        except OSError:
            pass
        try:
            os.remove("foo1.txt")
            os.remove("foo2.txt")
        except OSError:
            pass
         
    def test_loopCommandFull(self):
        """Use rename (mv) as a simple arbitrary command with testable result."""
        fu.do_foreach(self.files, ["mv", "***full***", "***full***.txt"], comments=False, progress=False)
        self.assertTrue(os.path.exists(self.files[0] + ".txt"))
        self.assertTrue(os.path.exists(self.files[1] + ".txt"))
        self.assertTrue(os.path.exists(self.files[2] + ".txt"))
            
    def test_loopCommandBase(self):
        """Use rename (mv) as a simple arbitrary command with testable result."""
        fu.do_foreach(self.files, ["mv", "./***file***", "./***file***.txt"], comments=False, progress=False)
        self.assertTrue(os.path.exists(self.files[0] + ".txt"))
        self.assertTrue(os.path.exists(self.files[1] + ".txt"))
        self.assertTrue(os.path.exists(self.files[2] + ".txt"))
         
    def test_loopCommandAlias(self):
        """Use rename (mv) as a simple arbitrary command with testable result."""
        fu.do_foreach(self.files, ["mv", "***full***", "./***alias***.txt"], comments=False, progress=False)
        self.assertTrue(os.path.exists(self.files[0] + ".txt"))
        self.assertTrue(os.path.exists(self.files[1] + ".txt"))
        self.assertTrue(os.path.exists(self.files[1] + ".txt"))
         
    def test_loopMainFiles(self):
        """Test command line use on files."""
        with open(os.devnull, 'w') as FNULL:
            command = ["python", script,
                       "T", self.files[0], self.files[1], 
                       "-C", "--loop", "S", "mv", "***full***", "***full***.txt"]
            subprocess.call(command, stdout=FNULL,  shell=False)
        self.assertTrue(os.path.exists(self.files[0] + ".txt"))
        self.assertTrue(os.path.exists(self.files[1] + ".txt"))
                  
    def test_loopMainFilesMixedSubstitution(self):
        """Test command line use on files."""
        with open(os.devnull, 'w') as FNULL:
            command = ["python", script,
                       "T", self.files[0], self.files[1], 
                       "-C", "--loop", "S", "mv", "***full***", "./***alias***.txt"]
            subprocess.call(command, stdout=FNULL,  shell=False)
        self.assertTrue(os.path.exists(self.files[0] + ".txt"))
        self.assertTrue(os.path.exists(self.files[1] + ".txt"))
      
    def test_loopMainRange(self):
        """Test command line use on a number range."""
        with open(os.devnull, 'w') as FNULL:
            command = ["python", script,
                       "T", "1:2", 
                       "-C", "--loop", "R", "mv", "testdofile***alias***", "testdofile***alias***.txt"]
            subprocess.call(command, stdout=FNULL,  shell=False)
        self.assertTrue(os.path.exists(self.files[0] + ".txt"))
        self.assertTrue(os.path.exists(self.files[1] + ".txt"))
          
    def test_loopMainToFiles(self):
        """Test main for outputting iterations to separate files."""
        with open(os.devnull, 'w') as FNULL:
            command = ["python", script,
                       "T", "1:2", 
                       "-C", "--loop", "R", "cat", "testdofile***alias***",
                       "-O", "./", "foo", ".txt"]
            subprocess.call(command, stdout=FNULL,  shell=False)
        self.assertTrue(os.path.exists("foo1.txt"))
        self.assertTrue(os.path.exists("foo2.txt"))
            
            
class Test_FilesEmpty(unittest.TestCase):
    def setUp(self):
        # Prepare two test files.
        self.empty = "emptytestfile"
        with open(self.empty,'w'):
            pass        
        self.notempty = "notemptytestfile"
        with open(self.notempty,'w') as f:
            f.write("foo bar")
               
    def tearDown(self):
        os.remove(self.empty)
        os.remove(self.notempty)
               
    def test_filesAreEmpty(self):
        """Test if empty and not empty files are identified correctly."""
        # Invalid as empty (Default)
        self.assertEqual(fu.are_empty( fu.FilesList([self.empty, self.notempty, _invalidPath]) ), 
                         fu.FilesList([self.empty, _invalidPath]))
        # Invalid as not empty
        self.assertEqual(fu.are_empty(fu.FilesList([self.empty, self.notempty, _invalidPath]), False), 
                         fu.FilesList([self.empty]))
             
    def test_filesEmptyMain(self):
        """Test command line use."""
        cwd = os.getcwd()
        command = ["python", script,
                    "T", self.empty, self.notempty,
                    "-C", "--probe", "are_empty"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read().rstrip()
        self.assertEqual(os.path.join(cwd, self.empty) +"\t"+ self.empty, result)
              
           
class Test_FilesMissing(unittest.TestCase):
    def setUp(self):
        # Create a file
        self.exists = "existingtestfile"
        with open(self.exists, 'w'):
            pass
                
    def tearDown(self):
        os.remove(self.exists)
                
    def test_filesDontExist(self):
        """ Test that valid and invalid paths are recognised correctly."""
        self.assertEqual(fu.dont_exist( fu.FilesList([self.exists, _invalidPath]) ), 
                         fu.FilesList([_invalidPath]))
                 
    def test_filesMissingMain(self):
        """Test command line use."""
        cwd = os.getcwd()
        command = ["python", script,
                    "T", self.exists, _invalidPath,
                    "-C", "--probe", "dont_exist"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read().rstrip()
        self.assertEqual(os.path.join(cwd, _invalidPath) +"\t"+ os.path.basename(_invalidPath), result)
           
           
class Test_FilesUnreadable(unittest.TestCase):
    def setUp(self):
        self.readable = "accesstestfile"
        with open(self.readable, 'w'):
            pass
        os.chmod(self.readable, 0400)
        self.unreadable = "inaccesstestfile"
        with open(self.unreadable, 'w'):
            pass
        os.chmod(self.unreadable, 0300)
                    
    def tearDown(self):
        os.chmod(self.readable, 0700)
        os.chmod(self.unreadable, 0700)
        os.remove(self.readable)
        os.remove(self.unreadable)
               
    def test_filesAreUnreadable(self):
        """Test read access to files."""
        self.assertEqual(fu.arent_readable( fu.FilesList([self.readable, self.unreadable, _invalidPath]) ), 
                         fu.FilesList([self.unreadable, _invalidPath]))
                  
    def test_filesUnreadableMain(self):
        """Test command line use."""
        cwd = os.getcwd()
        command = ["python", script,
                    "T", self.readable, self.unreadable,
                    "-C", "--probe", "arent_access"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read().rstrip()
        self.assertEqual(os.path.join(cwd, self.unreadable) +"\t"+ self.unreadable, result)
            
           
class Test_StringText(unittest.TestCase):
    def setUp(self):
        self.s = r"0123456789qwertyuiop[]asdfghjkl;'\\ zxcvbnm,./`!@$%^&*()_+}{POIUYTREWQASDFGHJKL:\"|?><MNBVCXZ~}"
              
    def test_stringIsText(self):
        """Test that plain text string evaluates returns True."""
        self.assertTrue(fu.istext(self.s))
              
    def test_stringIsntText(self):
        """Test that compressed text string returns False."""
        self.assertFalse(fu.istext(zlib.compress(self.s)))
           
        
class Test_FilesText(unittest.TestCase):
    def setUp(self):
        self.s = r"0123456789qwertyuiop[]asdfghjkl;'\\ zxcvbnm,./`!@$%^&*()_+}{POIUYTREWQASDFGHJKL:\"|?><MNBVCXZ~}"
        self.text = "texttestfile"
        with open(self.text, 'w') as f:
            f.write(self.s)
        self.empty = "emptytestfile"
        with open(self.empty, 'w'):
            pass
        self.binary = "binarytestfile"
        with open(self.binary, 'w') as f:
            f.write(zlib.compress(self.s))
                   
    def tearDown(self):
        os.remove(self.text)
        os.remove(self.empty)
        os.remove(self.binary)
           
    def test_filesArentText(self):
        """Test that files are recognised ocrrectly as text or binary."""
        self.assertEqual(fu.arent_text( fu.FilesList([self.text, self.binary, self.empty, _invalidPath]) ), 
                                       fu.FilesList([self.binary, _invalidPath]))
                     
    def test_main(self):
        """Test command line use."""
        cwd = os.getcwd()
        command = ["python", script,
                    "T", self.text, self.binary,
                    "-C", "--probe", "arent_text"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read().rstrip()
        self.assertEqual(os.path.join(cwd, self.binary) +"\t"+ self.binary, result)
            
           
class Test_SwapSubstrings(unittest.TestCase):
    def setUp(self):
        self.commastr= "11,12,13,14,15\n21,22,23,24,25\n31,32,33,34,35\n"
        self.mixstr= "11;12,13;14\t15\n21,22,23;24,25\n31\t32,33;34,35\n"
        self.tabstr = "11\t12\t13\t14\t15\n21\t22\t23\t24\t25\n31\t32\t33\t34\t35\n"
                     
    def test_swapDefaults(self):
        """Test default replacement of column separators."""
        self.assertEqual(fu.swap_substr([self.commastr]),
                          [self.tabstr])
                 
    def test_swapCustom(self):
        """Test custom replacement of column separators."""
        self.assertEqual(fu.swap_substr([self.tabstr], insep=["\t"], outsep=","),
                         [self.commastr])
                 
    def test_swapMixed(self):
        """Test replacing multiple string delimiters with a single one."""
        self.assertEqual(fu.swap_substr([self.mixstr], insep=[",",";"]),
                         [self.tabstr])
           
    def test_swapMultiple(self):
        """Test replacement of regex delimiter in multiple strings."""
        self.assertEqual(fu.swap_substr([self.commastr, self.mixstr], insep=[",|;"]),
                         [self.tabstr, self.tabstr])
              
      
class Test_SwapTabsInFiles(unittest.TestCase):
    def setUp(self):
        self.commastr= "11,12,13,14,15\n21,22,23,24,25\n31,32,33,34,35\n"
        self.csv = "csvtestfile"
        with open(self.csv, 'w') as f:
            f.write(self.commastr)
        self.mixstr= "11;12,13;14\t15\n21,22,23;24,25\n31\t32,33;34,35\n"
        self.msv = "msvtestfile"
        with open(self.msv, 'w') as f:
            f.write(self.mixstr)
        self.tabstr = "11\t12\t13\t14\t15\n21\t22\t23\t24\t25\n31\t32\t33\t34\t35\n"
        self.tsv = "tsvtestfile"
        with open(self.tsv, 'w') as f:
            f.write(self.tabstr)
            
    def tearDown(self):
        os.remove(self.csv)
        os.remove(self.tsv)
        os.remove(self.msv)
        try:
            os.remove("foo.csvtestfile.bar")
            os.remove("foo.msvtestfile.bar")
        except OSError:
            pass
                      
    def test_filesSwapDefaults(self):
        """Test default replacement of column separators."""
        self.assertEqual(fu.swap_strFiles([self.csv]),
                         [self.tabstr])
                  
    def test_filesSwapCustom(self):
        """Test custom replacement of column separators."""
        self.assertEqual(fu.swap_strFiles([self.tsv], insep=["\t"], outsep=","),
                         [self.commastr])
                  
    def test_filesSwapMixed(self):
        """Test replacing multiple string delimiters with a single one."""
        self.assertEqual(fu.swap_strFiles([self.msv], insep=[",",";"]),
                         [self.tabstr])
            
    def test_filesSwapMultiple(self):
        """Test replacement of regex delimiter in multiple files"""
        self.assertEqual(fu.swap_strFiles([self.csv, self.msv], insep=[",|;"]),
                         [self.tabstr, self.tabstr])
        
    def test_filesSwapMain(self):
        """Test command line use of swapping delimiters."""
        command = ["python", script,
                    "T", self.csv, self.msv,
                    "-C", "--swap", "\t",
                    "-s", ",", ";", ]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        self.assertEqual(result, self.tabstr + self.tabstr)
         
    def test_filesSwapMainToFiles(self):
        """Test command line use of swapping delimiters with individual outputs."""
        command = ["python", script,
                    "T", self.csv, self.msv,
                    "-C", "--swap", "\t",
                    "-s", ",", ";", 
                    "-O", "./", "foo.", ".bar"]
        subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        with open("foo.csvtestfile.bar") as f:
            self.assertEqual(f.read(), self.tabstr)
        with open("foo.msvtestfile.bar") as f:
            self.assertEqual(f.read(), self.tabstr)
          
           
class Test_CountFields(unittest.TestCase):
    def setUp(self):
        # Create minimal input files.
        self.vals = ["a","b","c","d","e"]
        self.comm = "# Comment."
        self.csv1 = "countcolstestfile1.csv"
        self.csv2 = "countcolstestfile2.csv"
        self.tsv = "countcolstestfile.tsv"
        with open(self.csv1, 'w') as f:
            f.write(self.comm + "\n")
            f.write(",".join(self.vals) + "\n")
        with open(self.csv2, 'w') as f:
            f.write(",".join(self.vals) + "," + ",".join(self.vals) + "\n")
            f.write(self.comm + "\n")
        with open(self.tsv, 'w') as f:
            f.write(self.comm + "\n")
            f.write("\t".join(self.vals) + "\n")
                
    def tearDown(self):
        os.remove(self.csv1)
        os.remove(self.csv2)
        os.remove(self.tsv)
        
    def test_countFieldsDefaults(self):
        """Test default behaviour."""
        self.assertEqual(fu.count_columns(fu.FilesList([self.tsv])), 
                         [len(self.vals)])
                
    def test_countFieldsCustom(self):
        """Test custom behaviour."""
        self.assertEqual(fu.count_columns(fu.FilesList([self.csv1, self.csv2]), colSep=[","]), 
                         [len(self.vals), 2* len(self.vals)])
                
    def test_countFieldsMain(self):
        """Test command line use."""
        command = ["python", script,
                   "T", self.csv1, self.tsv,
                   "-C", "--cntcols", 
                   "-s", ",", "\t"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        result = re.split("\t|\n", result)
        self.assertEqual(result[0], "5")
        self.assertEqual(result[3], "5")
        self.assertEqual(result[2], os.path.abspath(self.csv1))
        self.assertEqual(result[5], os.path.abspath(self.tsv)) 
     
     
class Test_PrepareDataframes(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame([[11,12,13,14],[21,22,23,24],[31,32,33,34]])
        self.ref = pd.DataFrame([[11,12,13,14],[21,22,23,24],[31,32,33,34]])
           
    def tearDown(self):
        pass
           
    def test_prepareDfDefaults(self):
        """Rename columns, don't drop first row, row labels"""
        df = fu.prepare_df(self.df)
        self.ref.columns = ["_"+ str(i) for i in range(0, self.ref.shape[1])]
        pd.util.testing.assert_frame_equal(df, self.ref)
              
    def test_prepareDfCustom(self):
        """Overriden column rename, drop first row, row labels, drop or don't drop index."""
        df = fu.prepare_df(self.df, myalias="foo", keyCol=1, keyhead="key", header=True, cols=[66,67,68,69])
        self.ref.columns = ["foo_66", "foo_67", "foo_68", "foo_69"]
        self.ref.set_index("foo_67", drop=False, inplace=True)
        self.ref.drop(self.ref.index.values.tolist()[0], axis=0, inplace=True)
        self.ref.index = [22,32]
        self.ref.index.name = "key"
        pd.util.testing.assert_frame_equal(df, self.ref)
             
           
class Test_GetColumnsFromFiles(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame([["11","12","13","14"],["21","22","23","24"],["31","32","33","34"]])
        self.df2 = pd.DataFrame([["11","12","13","14"],["31","32","33","34"],["21","22","23","24"]])
        self.df.index.name = None
        self.csv = "csvtestfile"
        self.df.to_csv(self.csv, header=False, index=False)
        self.csv2 = "csvtestfile2"
        self.df.to_csv(self.csv2, header=False, index=False)
        self.tsv = "tsvtestfile"
        self.df.to_csv(self.tsv, header=False, index=False, sep="\t")
        self.tsv2 = "tsvtestfile2"
        self.df2.to_csv(self.tsv2, header=False, index=False, sep="\t")
                      
    def tearDown(self):
        os.remove(self.csv)
        os.remove(self.csv2)
        os.remove(self.tsv)
        os.remove(self.tsv2)
        try:
            os.remove("foo.csvtestfile.bar")
            os.remove("foo.tsvtestfile.bar")
        except OSError:
            pass  
                      
    def test_getColumnsDefaults(self):
        """Test getting defaults from one file (first column, no header, '\\t' separator)."""
        # The reference value must be a dataframe with columns labeled with the originating file name and column.
        slice = self.df.iloc[:,0:1] # Without the slice it would be a Series object, not DataFrame.
        slice.columns = [self.tsv]
        result = fu.get_columns([self.tsv])[0]
        pd.util.testing.assert_frame_equal(result, slice)
                 
    def test_getColumnsCustom(self):
        """Test getting a custom slice from multiple files."""
        # Set up reference values with correct labels and correct index.
        idxs = [1,2,3,1,2,3]
        slice = self.df.iloc[1:,idxs]
        slice.columns = [self.csv +'_'+str(idxs[0]), self.csv +'_'+str(idxs[1]), self.csv +'_'+str(idxs[2]),
                         self.csv2 +'_'+str(idxs[3]), self.csv2 +'_'+str(idxs[4]), self.csv2 +'_'+str(idxs[5])]
        result = fu.get_columns([self.csv, self.csv2], cols=range(1,4), colSep=[","], header=True)[0]
        pd.util.testing.assert_frame_equal(result,slice)
        result = fu.get_columns([self.csv, self.csv2], cols=["1:3"], colSep=[","], header=True)[0]
        pd.util.testing.assert_frame_equal(result,slice)
           
    def test_getColumnsIndexed(self):
        """Test that row indexes are handled correctly."""
        flist = fu.FilesList([self.csv, self.tsv2])
        # Set up reference values with correct labels and correct index.
        idxs = [3,2,3,2]
        slice = self.df.iloc[1:,idxs]
        slice.columns = [self.csv +'_'+str(idxs[0]), self.csv +'_'+str(idxs[1]), self.tsv2 +'_'+str(idxs[2]), self.tsv2 +'_'+str(idxs[3])]
        slice.index = self.df.iloc[1:,0]
        slice.index.name = str(self.df.iloc[0,0])
        result = fu.get_columns(flist, colSep=[",", "\t"], cols=[3,2], index=0, header=True)[0]
        pd.util.testing.assert_frame_equal(result,slice)
                 
    def test_getColumnsNoMerging(self):
        """Test that results can be returned un-merged."""
        flist = fu.FilesList([self.csv, self.csv2])
        result = fu.get_columns(flist, colSep=[","], cols=[3,2], index=0, header=True, merge=False)
        # Set up reference values with correct labels and correct index.
        idxs = [3,2]
        slice = self.df.iloc[1:,idxs]
        slice.columns = [self.csv +'_'+str(idxs[0]), self.csv +'_'+str(idxs[1])]
        slice.index = self.df.iloc[1:,0]
        slice.index.name = str(self.df.iloc[0,0])
        pd.util.testing.assert_frame_equal(result[0], slice)
        # Set up next reference.
        slice.columns = [self.csv2 +'_'+str(idxs[0]), self.csv2 +'_'+str(idxs[1])]
        slice.index = self.df.iloc[1:,0]
        slice.index.name = str(self.df.iloc[0,0])
        pd.util.testing.assert_frame_equal(result[1], slice)
               
    def test_getColumnsMain(self):
        """Test command line use."""
        command = ["python", script,
                   "T", self.csv, self.tsv,
                   "-C", "--cols", "3", "2", 
                   "-s", ",", "\t"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        self.assertEqual(result.rstrip("\n"), 
                         "csvtestfile_3,csvtestfile_2,tsvtestfile_3,tsvtestfile_2\n14,13,14,13\n24,23,24,23\n34,33,34,33")
        # And with mixed format
        command = ["python", script,
                   "T", self.csv, self.tsv,
                   "-C", "--cols", "3,0", "1:2", 
                   "-s", ",", "\t"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        self.assertEqual(result.rstrip("\n"), 
                         "csvtestfile_3,csvtestfile_0,csvtestfile_1,csvtestfile_2,\
                          tsvtestfile_3,tsvtestfile_0,tsvtestfile_1,tsvtestfile_2,\
                          \n14,11,12,13,14,11,12,13\n24,21,22,23,24,21,22,23\n34,31,32,33,34,31,31,33")
        # And one without relabel and headers.
        command = ["python", script,
                   "T", self.csv, self.tsv,
                   "-C", "--cols", "3", "2", 
                   "-s", ",", "\t", "-rl"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        self.assertEqual(result.rstrip("\n"), 
                         "24,23,24,23\n34,33,34,33")
        # And one with index
        command = ["python", script,
                   "T", self.csv, self.tsv2,
                   "-C", "--cols", "3", "2", 
                   "-s", ",", "\t", "-rl", "-i"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        self.assertEqual(result.rstrip("\n"), 
                         "21,24,23,24,23\n31,34,33,34,33")
                  
    def test_getColumnsMainToFiles(self):
        """Test command line use with output to individual files."""
        command = ["python", script,
                   "T", self.csv, self.tsv,
                   "-C", "--cols", "3", "2", 
                   "-s", "\t", ",",
                   "-O", "./", "foo.", ".bar"]
        subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        with open("foo.csvtestfile.bar") as f:
            self.assertEqual(f.read().rstrip("\n"), "csvtestfile_3\tcsvtestfile_2\n14\t13\n24\t23\n34\t33")
        with open("foo.tsvtestfile.bar") as f:
            self.assertEqual(f.read().rstrip("\n"), "tsvtestfile_3\ttsvtestfile_2\n14\t13\n24\t23\n34\t33")
             
            
class Test_GetColumnsFromFileManually(unittest.TestCase):
    def setUp(self):
        self.good=pd.DataFrame([["1","2","3"],["4","5","6"],["7","8","9"]])
        self.bad="1,2\n4\n7,8,9,10\n"
        self.badref=pd.DataFrame([["1","2","IDXERROR","IDXERROR"],["4","IDXERROR","IDXERROR","IDXERROR"],["7","8","9","10"]])
        self.csvg="csvtestfileG"
        self.csvb="csvtestfileB"
        self.tsvg="tsvtestfileG"
        self.tsvb="tsvtestfileB"
        self.good.to_csv(self.tsvg, sep="\t", header=False, index=False)
        self.good.to_csv(self.csvg, header=False, index=False)
        with open(self.csvb, 'w') as f:
            f.write(self.bad)
        with open(self.tsvb, 'w') as f:
            f.write(self.bad.replace(",", "\t"))
                    
    def tearDown(self):
        os.remove(self.csvb)
        os.remove(self.csvg)
        os.remove(self.tsvb)
        os.remove(self.tsvg)
                    
    def test_getColumnsManuallyDefaults(self):
        """Test default behaviour."""
        self.good.columns = [ self.tsvg for i in [0,1,2]]
        result = fu.get_columns_manual(self.tsvg)
        slice = self.good.iloc[:,0:1]
        pd.util.testing.assert_frame_equal(result,slice)
        self.badref.columns = [ self.tsvb for i in [0,1,2,3]]
        slice = self.badref.iloc[:,0:1]
        result = fu.get_columns_manual(self.tsvb)
        pd.util.testing.assert_frame_equal(result, slice)

    def test_getColumnsManuallyCustom(self):
        """Test parametrized behaviour."""
        self.good.columns = [ self.csvg +"_"+ str(i) for i in [0,1,2]]
        pd.util.testing.assert_frame_equal(fu.get_columns_manual(self.csvg, cols=[1,2], colSep=[","]), 
                                           self.good.iloc[:,[1,2]])
        pd.util.testing.assert_frame_equal(fu.get_columns_manual(self.csvg, cols=["1:2"], colSep=[","]), 
                                           self.good.iloc[:,[1,2]])
        self.badref.columns = [ self.csvb +"_"+ str(i) for i in [0,1,2,3]]
        pd.util.testing.assert_frame_equal(fu.get_columns_manual(self.csvb, cols=[1,2], colSep=[","]), 
                                           self.badref.iloc[:,[1,2]])
        self.badref.columns = [ self.csvb for i in [0,1,2,3]]
        result = fu.get_columns_manual(self.csvb, cols=[3], colSep=[","])
        slice = self.badref.iloc[:,3:4]
        pd.util.testing.assert_frame_equal(result, slice)
#         
    def test_getColumnsIndexed(self):
        """Test that row indexes are handled correctly."""
        result = fu.get_columns_manual(self.tsvg, index=0, cols=[2,1])
        labels = [ self.tsvg +"_0", self.tsvg +"_1", self.tsvg +"_2"]
        self.good.columns = labels
        self.good.set_index(labels[0], inplace=True, drop=False)
        self.good.index.name = "row_ID"
        slice = self.good.loc[:,[labels[2], labels[1]]]
        pd.util.testing.assert_frame_equal(result, slice )
        # What if the index is among the columns i want?
        result = fu.get_columns_manual(self.tsvg, index=0, cols=[0,1])
        slice = self.good.loc[:,[labels[0], labels[1]]]
        pd.util.testing.assert_frame_equal(result, slice)
               
       
class Test_GetRandomColumnsFromFiles(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame([[11,12,13,14],[21,22,23,24],[31,32,33,34]])
        self.csv1 = "csvtestfile1"
        self.df.to_csv(self.csv1, header=False, index=False)
        self.csv2 = "csvtestfile2"
        self.df.to_csv(self.csv2, header=False, index=False)
        self.tsv = "tsvtestfile"
        self.df.to_csv(self.tsv, header=False, index=False, sep="\t")
        self.df2 = pd.DataFrame([[11,12,13,14],[31,32,33,34],[21,22,23,24]])
        self.tsv2 = "tsvtestfile2"
        self.df2.to_csv(self.tsv2, header=False, index=False, sep="\t")
         
    def tearDown(self):
        os.remove(self.csv1)
        os.remove(self.csv2)
        os.remove(self.tsv)
        os.remove(self.tsv2)
        try:
            os.remove("foo.csvtestfile1.bar")
            os.remove("foo.tsvtestfile.bar")
        except OSError:
            pass  
                 
    def test_getRandColumnsDefaults(self):
        """Test getting defaults from one file (one column, no header, '\\t' separator)."""
        # There is no way to know in advance which column will be returned each time..
        # I can only check that the correct number of columns is returned.
        df = fu.get_random_columns([self.tsv])[0]
        self.assertEqual(len(df.columns.values.tolist()), 1)
                  
    def test_getRandColumnsCustom(self):
        """Test getting a custom slice from multiple files."""
        # Again, can only test the correct number of columns is returned.
        flist = fu.FilesList([self.csv1, self.csv2])
        df = fu.get_random_columns(flist, colSep=[","], k=3)[0]
        self.assertEqual(len(df.columns.values.tolist()), 3 * len(flist))
                  
    def test_getRandColumnsWithRowIndex(self):
        """Test that row indexes are handled correctly."""
        # Again, can only test the correct number of columns is returned.
        flist = fu.FilesList([self.csv1, self.csv2])
        df = fu.get_random_columns(flist, colSep=[","], k=2, index=0, header=True)[0]
        self.df.set_index(self.df.index.values.tolist()[0], drop=False, inplace=True)
        self.df.drop(self.df.index.values.tolist()[0], axis=0, inplace=True)
        self.df.index.name = "row_ID"
        # Check number of columns.
        self.assertEqual(len(df.columns.values.tolist()), 2 * len(flist))
        # Check the index.
        self.assertEqual(df.index.values.tolist(), self.df.index.values.tolist())
        # Check that the index header appears only once.
        self.assertTrue(df.index.name not in df.columns.values.tolist())
                
    def test_getRandomColumnsNoMerging(self):
        """Test that it can return individual results."""
        # Again, can only test the correct number of columns is returned.
        result = fu.get_random_columns([self.csv1, self.tsv], colSep=["\t",","], 
                                       k=2, merge=False)
        # Check number of columns.
        self.assertEqual(len(result), 2)
        self.assertEqual(len(result[0].columns.values.tolist()), 2)
        self.assertEqual(len(result[1].columns.values.tolist()), 2)
        # Also try with index.
        result = fu.get_random_columns([self.csv1, self.tsv], colSep=["\t",","], 
                                       k=2, merge=False,  index=0)
        self.df.set_index(self.df.index.values.tolist()[0], drop=False, inplace=True)
        self.df.index.name = "row_ID"
        # Check number of columns.
        self.assertEqual(len(result), 2)
        self.assertEqual(len(result[0].columns.values.tolist()), 2)
        self.assertEqual(len(result[1].columns.values.tolist()), 2)
        # Check the index.
        self.assertEqual(result[0].index.values.tolist(), self.df.index.values.tolist())
        # Check that the index header appears only once.
        self.assertTrue(result[0].index.name not in result[0].columns.values.tolist())
        self.assertTrue(result[1].index.name not in result[1].columns.values.tolist())
               
    def test_getRandomColumnsMain(self):
        """Test command line use."""
        command = ["python", script,
                   "T", self.csv1, self.tsv,
                   "-C", "--rndcols", "2", 
                   "-s", ",", "\t"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        lines = result.rstrip("\n").split("\n")
        self.assertEqual(len(lines), 4)
        fields = lines[0].split(",")
        self.assertEqual(len(fields), 4) 
        # And one without relabel and headers.
        command = ["python", script,
                   "T", self.csv1, self.tsv,
                   "-C", "--rndcols", "2", 
                   "-s", ",", "\t", "-rl"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        lines = result.rstrip("\n").split("\n")
        self.assertEqual(len(lines), 2)
        fields = lines[0].split(",")
        self.assertEqual(len(fields), 4) 
        # And one with index
        command = ["python", script,
                   "T", self.csv1, self.tsv2,
                   "-C", "--rndcols", "2", 
                   "-s", ",", "\t", "-rl", "-i"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        lines = result.rstrip("\n").split("\n")
        self.assertEqual(len(lines), 2)
        fields = lines[0].split(",")
        self.assertEqual(len(fields), 5) 
        self.assertEqual(fields[0], "21")
                 
              
    def test_getRandomColumnsMainToFiles(self):
        """Test command line use with output to individual files."""
        # I already know by now if the index and header and relabeling works.
        # I just have to test that it can go to files.
        # I can again only test the structure, not the actual values.
        command = ["python", script,
                   "T", self.csv1, self.tsv,
                   "-C", "--rndcols", "2", 
                   "-s", "\t", ",",
                   "-O", "./", "foo.", ".bar"]
        subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        with open("foo.csvtestfile1.bar") as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 4)
            fields = lines[0].split("\t")
            self.assertEqual(len(fields), 2)
        with open("foo.tsvtestfile.bar") as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 4)
            fields = lines[0].split("\t")
            self.assertEqual(len(fields), 2)
    
              
class Test_appendColumns(unittest.TestCase):
    def setUp(self):
        self.csv1 = "csvtestfile1"
        self.csv2 = "csvtestfile2"
        self.tsv = "tsvtestfile"
        self.df = pd.DataFrame([["a","11","12","13","14"],
                                ["b","21","22","23","24"],
                                ["c","31","32","33","34"]])
        self.df.to_csv(self.tsv, header=False, index=False, sep="\t")
        self.df1 = pd.DataFrame([["Z","A","B","C","D"],
                                 ["a","11","12","13","14"],
                                 ["b","21","22","23","24"],
                                 ["c","31","32","33","34"]])
        self.df1.to_csv(self.csv1, header=False, index=False)
        self.df2 = pd.DataFrame([["z","e","f"],
                                 ["c","35","36"],
                                 ["a","15","16"],
                                 ["b","25","26"]])
        self.df2.to_csv(self.csv2, header=False, index=False)
           
                 
    def tearDown(self):
        os.remove(self.csv1)
        os.remove(self.csv2)
        os.remove(self.tsv)
                 
    def test_appendDefaults(self):
        """Test key-oblivious appending of columns using default parameters."""
        self.df.columns = [self.tsv +"_"+ str(i) for i in self.df.columns.values.tolist()]
        df = fu.append_columns([self.tsv, self.tsv])
        pd.util.testing.assert_frame_equal(df.iloc[:,0:5], self.df)
        pd.util.testing.assert_frame_equal(df.iloc[:,5:10], self.df)
            
    def test_appendDefaultsIndexed(self):
        """Test key-aware appending of columns using default parameters."""
        ref = pd.DataFrame([["a","11","12","13","14","11","12","13","14"],
                            ["b","21","22","23","24","21","22","23","24"],
                            ["c","31","32","33","34","31","32","33","34"]])
        ref.columns = ["row_ID","tsvtestfile_1","tsvtestfile_2","tsvtestfile_3","tsvtestfile_4",
                       "tsvtestfile_1","tsvtestfile_2","tsvtestfile_3","tsvtestfile_4"]
        ref.set_index("row_ID", drop=True, inplace=True)
        df = fu.append_columns([self.tsv, self.tsv], index=0)
        pd.util.testing.assert_frame_equal(df, ref)
                      
    def test_appendCustom(self):
        """Test key-oblivious appending of columns using custom parameters."""
        ref = pd.DataFrame([["Z","A","B","C","D","z","e","f"],
                            ["a","11","12","13","14","c","35","36"],
                            ["b","21","22","23","24","a","15","16"],
                            ["c","31","32","33","34","b","25","26"]])
        ref.drop(0, axis=0, inplace=True)
        ref.columns = ["csvtestfile1_0","csvtestfile1_1","csvtestfile1_2","csvtestfile1_3",
                       "csvtestfile1_4","csvtestfile2_0","csvtestfile2_1","csvtestfile2_2"]
        df = fu.append_columns([self.csv1, self.csv2], colSep=[","], header=True)
        pd.util.testing.assert_frame_equal(df, ref)
            
    def test_appendCustomIndexed(self):
        """Test key-aware appending of columns using custom parameters."""
        ref = pd.DataFrame([["a","11","12","13","14","15","16"],
                            ["b","21","22","23","24","25","26"],
                            ["c","31","32","33","34","35","36"]])
        ref.columns = ["Z","csvtestfile1_1","csvtestfile1_2","csvtestfile1_3",
                       "csvtestfile1_4","csvtestfile2_1","csvtestfile2_2"]
        ref.set_index("Z", drop=True, inplace=True)
        df = fu.append_columns([self.csv1, self.csv2], colSep=[","], header=True, index=0)
        pd.util.testing.assert_frame_equal(df, ref)
                     
    def test_appendMain(self):
        """Test the function and syntax of using through main."""
        command = ["python", script,
                   "T", self.csv1, self.csv2,
                   "-C", "--appnd",
                   "-s", ","]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        ref = "csvtestfile1_0,csvtestfile1_1,csvtestfile1_2,csvtestfile1_3,csvtestfile1_4,csvtestfile2_0,csvtestfile2_1,csvtestfile2_2\nZ,A,B,C,D,z,e,f\na,11,12,13,14,c,35,36\nb,21,22,23,24,a,15,16\nc,31,32,33,34,b,25,26"
        self.assertEqual(result.rstrip("\n"), ref)
        # Indexed, without labels or relabeling
        command = ["python", script,
                   "T", self.csv1, self.csv2,
                   "-C", "--appnd",
                   "-s", ",", "-i", "-rl"]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        ref = "a,11,12,13,14,15,16\nb,21,22,23,24,25,26\nc,31,32,33,34,35,36"
        self.assertEqual(result.rstrip("\n"), ref)
      
      
class Test_GetCrosspointsFromFiles(unittest.TestCase):
    def setUp(self):
        self.ref = pd.DataFrame([["Z","A","B","C","D"],
                            ["a","11","12","13","14"],
                            ["b","21","22","23","24"],
                            ["c","31","32","33","34"]])
        self.csv = "csvtestfile"
        self.tsv = "tsvtestfile"
        self.ref.to_csv(self.csv, header=False, index=False)
        self.ref.to_csv(self.tsv, sep="\t", header=False, index=False)
           
    def tearDown(self):
        os.remove(self.csv)
        os.remove(self.tsv)
           
    def test_crosspointsDefault(self):
        df = fu.get_crosspoints([self.tsv, self.tsv])[0]
        self.assertEqual(df.iloc[0,0], self.ref.iloc[0,0])
        self.assertEqual(df.iloc[0,1], self.ref.iloc[0,0])
        self.assertEqual(df.shape, (1,2))
       
    def test_crosspointsCustom(self):
        "Test getting the cells at intersections of selected rows and columns."
        results = fu.get_crosspoints([self.csv, self.csv], cols=[1,3], rows=[2,0], colSep=[","], header=True, index=3, merge=False)
        self.assertEqual(results[0].iloc[0,0], self.ref.iloc[3,1])
        self.assertEqual(results[1].iloc[0,0], self.ref.iloc[3,1])
        self.assertEqual(results[0].iloc[0,1], self.ref.iloc[3,3])
        self.assertEqual(results[1].iloc[0,1], self.ref.iloc[3,3])
        self.assertEqual(results[0].iloc[1,0], self.ref.iloc[1,1])
        self.assertEqual(results[1].iloc[1,0], self.ref.iloc[1,1])
        self.assertEqual(results[0].iloc[1,1], self.ref.iloc[1,3])
        self.assertEqual(results[1].iloc[1,1], self.ref.iloc[1,3])
        self.assertEqual(results[0].shape, (2,2))
        self.assertEqual(results[1].shape, (2,2))
          
 
class Test_GetValuesSet(unittest.TestCase):
    def setUp(self):
        self.ref = pd.DataFrame([["M","N","P","Q","Z"],
                            ["a","11","12","13","13"],
                            ["b","21","22","23","24"],
                            ["c","31","32","33","34"],
                            ["d","31","32","33","44"]])
        self.csv = os.path.abspath("csvtestfile")
        self.tsv = os.path.abspath("tsvtestfile")
        self.ref.to_csv(self.csv, header=False, index=False)
        self.ref.to_csv(self.tsv, sep="\t", header=False, index=False)
         
    def tearDown(self):
        os.remove(self.csv)
        os.remove(self.tsv)
         
    def test_getValSetDefaults(self):
        """Collection of non-redundant values with defaults."""
        nest = fu.get_valuesSet([self.tsv])
        self.assertEqual(nest[0], set(["M","N","P","Q","Z"]))
         
    def test_getValSetAxis(self):
        """Collection of non-redundant values in either axis."""
        nest = fu.get_valuesSet([self.tsv], colSep=["\t"], axis='c', index=2)
        self.assertEqual(nest[0], set(["P","12","22","32"]))
        nest = fu.get_valuesSet([self.tsv], colSep=["\t"], axis='r', index=0)
        self.assertEqual(nest[0], set(["M","N","P","Q","Z"]))
     
    def test_getValSetMultifile(self):
        """Collection of non-redundant values with multiple input."""
        nest = fu.get_valuesSet([self.tsv, self.csv], colSep=[",","\t"], axis='c', index=2)
        self.assertEqual(nest[0], set(["P","12","22","32"]))
        self.assertEqual(nest[0], nest[1])
     
    def test_getValSetMode(self):
        """Collection of non-redundant values with repetition filtering."""
        nest = fu.get_valuesSet([self.tsv], axis='c', index=3, filter='a')
        self.assertEqual(nest[0], set(["Q","13","23","33"]))
        nest = fu.get_valuesSet([self.tsv], axis='c', index=3, filter='u')
        self.assertEqual(nest[0], set(["Q","13","23"]))
        nest = fu.get_valuesSet([self.tsv], axis='c', index=3, filter='r')
        self.assertEqual(nest[0], set(["33"]))
 
    def test_getValSetErrors(self):
        """Collection of non-redundant values with invalid parameters."""
        self.assertRaises(ValueError, fu.get_valuesSet, [self.csv], axis='f')
        self.assertRaises(ValueError, fu.get_valuesSet, [self.csv], filter='f')
 
    def test_appendMain(self):
        """Collection of non-redundant values from Main."""
        command = ["python", script,
                   "T", self.csv,
                   "-C", "--valset", 'r', '1', 'a',
                   "-s", ","]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        self.assertEqual(result.rstrip("\n"), "\t".join([self.csv, str(set(["a","11","12","13"]))]))
        command = ["python", script,
                   "T", self.csv,
                   "-C", "--valset", 'c', '1', 'u',
                   "-s", ","]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        self.assertEqual(result.rstrip("\n"), "\t".join([self.csv, str(set(["N","11","21"]))]))
        command = ["python", script,
                   "T", self.csv,
                   "-C", "--valset", 'c', '1', 'r',
                   "-s", ","]
        result = subprocess.Popen(command, stdout=subprocess.PIPE,  shell=False).stdout.read()
        self.assertEqual(result.rstrip("\n"), "\t".join([self.csv, str(set(["31"]))]))




if __name__ == "__main__":
    unittest.main()


#EOF