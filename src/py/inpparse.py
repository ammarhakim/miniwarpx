r"""Provides classes to read and parse WarpX chapter-section style
input files.
"""

import sre
from optparse import OptionParser

class Section:
    r"""Represents a section of a WarpX input file."""
    
    def __init__(self, name):
        self.name = name # section name
        self.data = {} # key -> value

    def getValue(self, key, default=None):
        r"""getValue(key:str, default:type) -> str

        Returns value associated with `key` if it exists, or `default`
        otherwise.
        """
        if key in self.data:
            return self.data[key]
        else:
            return default # no such key exists

    def addValue(self, key, value):
        r"""addValue(key:str, value:str) -> None

        Adds `key` `value` pair to the section.
        """

        self.data[key] = value

    def varNames(self):
        r"""varNames() -> [] str

        Returns list of variables defined in section
        """
        return self.data.keys()

class Chapter:
    r"""Represents a chapter (i.e. the complete input file)."""
    
    def __init__(self):
        self.sections = {} # sectionName -> Section

    def getSection(self, sectionName):
        r"""getSection(sectionName:str) -> Section

        Returns Section with name `sectionName`.
        """
        if sectionName in self.sections:
            return self.sections[sectionName]
        else:
            return None # no such section exists

    def addSection(self, sectionName, section):
        r"""addSection(sectionName:str, section:Section) -> None

        Adds section with specified name into chapter.
        """
        self.sections[sectionName] = section

    def sectionNames(self):
        r"""sectionNames() -> [] str

        Returns list of section names in chapter.
        """
        return self.sections.keys()
        

class InpParser:

    # matches "[SectionName]"
    reSection  = sre.compile('\[(?P<sn>\S+)\]')
    # matches "key = value" 
    reKeyValue = sre.compile('(?P<key>\w+)[ \t]*=[ \t]*(?P<value>\S+)')
    
    def __init__(self, fname):
        self.fname = fname # input file name

    def parse(self):
        r"""parse() -> Chapter

        Parses input file and returns chapter object for the complete
        input file."""

        chapter = Chapter()

        # open input file for reading
        fp = open(self.fname, 'r')
        data = fp.read().split('\n') # split into lines

        currSection = None
        for d in data:
            # check if line is section header
            mo = self.reSection.search(d)
            if mo:
                # yes, it is a section header: extract name and add it to chapter
                sn = mo.group('sn')
                section = Section(sn)
                chapter.addSection(sn, section)
                currSection = section
            else:
                # no it is not: check if it is key=value pair
                mo = self.reKeyValue.search(d)
                if mo:
                    # yes, it is: extract key and pair
                    key = mo.group('key')
                    value = mo.group('value')
                    currSection.addValue(key,value)
                else:
                    # it is neither section nor data line so ignore it
                    pass

        return chapter

class InpCmdLineOptions:
    
    r"""Reads a chapter-section style template file constructing the list of
    variables that can be specified on the command prompt.
    """

    # matches "$var"
    varRe = sre.compile('\$(?P<var>\w+)')

    def __init__(self, fname):
        self.fname = fname # file name
        self.vars = {} # dictionary of variable -> value

    def parse(self):
        r"""parse() -> None
        """

        fp = open(self.fname, 'r')
        data = fp.read()
        values = self.varRe.findall(data) # findall all $vars in file

        # loop over values, inserting them in the dictionary
        for v in values:
            self.vars[v] = None

    def getVars(self):
        return self.vars.keys()
