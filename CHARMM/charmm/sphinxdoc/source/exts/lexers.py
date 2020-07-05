from pygments.lexer import RegexLexer, bygroups
from pygments.token import *
import re

# These determine how to tokenize CHARMM commands
ChCommand = Number
SelectKey = String
Selection = Name.Variable.Global
UserVarbl = Name.Variable.Class
ChrmVarbl = Name.Variable.Class
FileTitle = Comment.Special

class ChmLexer(RegexLexer):
    ''' Defines a lexer to enable syntax highlighting in CHARMM input files. '''
    name = 'CHARMM INPUT'
    aliases = ['chm']
    filenames = ['*.inp']

    flags = re.IGNORECASE | re.MULTILINE | re.DOTALL

    tokens = {
        'root': [
            # * marks a title; read until the end of this line
            (r'\*.*?$', FileTitle),        
            # Handle the ! comments
            (r'\s*\!.*?$', Comment),
            # The first word on a line is a command; hop into line subparser
            (r'\s*\w+', ChCommand, 'line'),
            # Everything else is just text
            (r'.', Text),
        ],
        'line': [
            # Handle the ! comments
            (r'\s*\!.*?$', Comment),
            # If we find a newline with hyphen (and possibly a comment) before it, it
            # is a line continuation, so we just grab the comment and carry on
            (r'(-)(\!.*?)?(\n)', bygroups(Text, Comment, Text)),
            # If we find a newline (no hyphen before it), with nothing else (but a comment)
            # at the end of it, flag it and go back up
            (r'(\s*\!.*?)?(\n)', bygroups(Comment, Text), '#pop'),
            # Look for atom selection; hop into the select parser if we find one
            (r'SELE\w*', SelectKey, 'select'),
            # CHARMM internal variables
            (r'\?\w+', ChrmVarbl),
            # CHARMM user variables
            (r'@\w+', UserVarbl),
            # Everything else is just text
            (r'.', Text),
        ],
        'select': [
            # If we find a newline with hyphen (and possibly a comment) before it, it
            # is a line continuation, so we just grab the comment and carry on
            (r'(-)(\!.*?)?(\n)', bygroups(Selection, Comment, Selection)),
            # If we find an end command, return
            (r'END', SelectKey, '#pop'),
            # Everything else is part of the selection
            (r'.', Selection),
        ]
    }


class CrdLexer(RegexLexer):
    ''' Defines a lexer to enable syntax highlighting in CHARMM coordinate files. '''
    name = 'CHARMM Coordinates'
    aliases = ['crd']
    filenames = ['*.crd']

    flags = re.IGNORECASE | re.MULTILINE | re.DOTALL

    tokens = {
        'root': [
            # * marks a title; read until the end of this line
            (r'\*.*?$', FileTitle),        
            # Handle the ! comments
            (r'\s*\!.*?$', Comment),
            # Everything else is just text
            (r'.', Text),
        ],
    }


class RTFLexer(RegexLexer):
    ''' Defines a lexer to enable syntax highlighting in CHARMM RTF files. '''
    name = 'CHARMM RTF'
    aliases = ['rtf']
    filenames = ['*.rtf']

    flags = re.IGNORECASE | re.MULTILINE | re.DOTALL

    tokens = {
        'root': [
            # * marks a title; read until the end of this line
            (r'\*.*?$', FileTitle),        
            # Handle the ! comments
            (r'\s*\!.*?$', Comment),
            # Everything else is just text
            (r'.', Text),
        ],
    }

class PRMLexer(RegexLexer):
    ''' Defines a lexer to enable syntax highlighting in CHARMM PRM files. '''
    name = 'CHARMM PRM'
    aliases = ['prm']
    filenames = ['*.prm']

    flags = re.IGNORECASE | re.MULTILINE | re.DOTALL

    tokens = {
        'root': [
            # * marks a title; read until the end of this line
            (r'\*.*?$', FileTitle),        
            # Handle the ! comments
            (r'\s*\!.*?$', Comment),
            # Everything else is just text
            (r'.', Text),
        ],
    }


class PSFLexer(RegexLexer):
    ''' Defines a lexer to enable syntax highlighting in CHARMM PSF files. '''
    name = 'CHARMM PSF'
    aliases = ['psf']
    filenames = ['*.psf']

    flags = re.IGNORECASE | re.MULTILINE | re.DOTALL

    tokens = {
        'root': [
            # * marks a title; read until the end of this line
            (r'\*.*?$', FileTitle),        
            # Handle the ! comments
            (r'\s*\!.*?$', Comment),
            # Everything else is just text
            (r'.', Text),
        ],
    }
