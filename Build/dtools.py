# #START_LICENSE###########################################################
#
# Parts of the code used here is modified from the Environment for Tree
# Exploration program (ETE).  http://etetoolkit.org
#
# Both this code and ETE is free software:
# you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ETE is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# <http://www.gnu.org/licenses/>.
#
#
# #END_LICENSE#############################################################

import copy, os
import argparse

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt

import random
import copy
import itertools
from collections import deque
from hashlib import md5
from functools import cmp_to_key

#import six
#from six.moves import (cPickle, map, range, zip)


import re

ITERABLE_TYPES = set([list, set, tuple, frozenset])

# Regular expressions used for reading newick format
_ILEGAL_NEWICK_CHARS = ":;(),\[\]\t\n\r="
_NON_PRINTABLE_CHARS_RE = "[\x00-\x1f]+"

_NHX_RE = "\[&&NHX:[^\]]*\]"
_FLOAT_RE = "\s*[+-]?\d+\.?\d*(?:[eE][-+]?\d+)?\s*"
# _FLOAT_RE = "[+-]?\d+\.?\d*"
# _NAME_RE = "[^():,;\[\]]+"
_NAME_RE = "[^():,;]+?"

# thanks to: http://stackoverflow.com/a/29452781/1006828
_QUOTED_TEXT_RE = r"""((?=["'])(?:"[^"\\]*(?:\\[\s\S][^"\\]*)*"|'[^'\\]*(?:\\[\s\S][^'\\]*)*'))"""
# _QUOTED_TEXT_RE = r"""["'](?:(?<=")[^"\\]*(?s:\\.[^"\\]*)*"|(?<=')[^'\\]*(?s:\\.[^'\\]*)*')""]"]"""
# _QUOTED_TEXT_RE = r"""(?=["'])(?:"[^"\\]*(?:\\[\s\S][^"\\]*)*"|'[^'\\]*(?:\\[\s\S][^'\\]*)*')]"]")"]"""

_QUOTED_TEXT_PREFIX = 'ete3_quotref_'

DEFAULT_DIST = 1.0
DEFAULT_NAME = ''
DEFAULT_SUPPORT = 1.0
FLOAT_FORMATTER = "%0.6g"
# DIST_FORMATTER = ":"+FLOAT_FORMATTER
NAME_FORMATTER = "%s"


def set_float_format(formatter):
    ''' Set the conversion format used to represent float distances and support
    values in the newick representation of trees.

    For example, use set_float_format('%0.32f') to specify 32 decimal numbers
    when exporting node distances and bootstrap values.

    Scientific notation (%e) or any other custom format is allowed. The
    formatter string should not contain any character that may break newick
    structure (i.e.: ":;,()")

    '''
    global FLOAT_FORMATTER
    FLOAT_FORMATTER = formatter
    # DIST_FORMATTER = ":"+FLOAT_FORMATTER


# Allowed formats. This table is used to read and write newick using
# different convenctions. You can also add your own formats in an easy way.
#
#
# FORMAT: [[LeafAttr1, LeafAttr1Type, Strict?], [LeafAttr2, LeafAttr2Type, Strict?],\
#    [InternalAttr1, InternalAttr1Type, Strict?], [InternalAttr2, InternalAttr2Type, Strict?]]
#
# Attributes are placed in the newick as follows:
#
# .... ,LeafAttr1:LeafAttr2)InternalAttr1:InternalAttr2 ...
#
#
#           /-A
# -NoName--|
#          |          /-B
#           \C-------|
#                    |          /-D
#                     \E-------|
#                               \-G
#
# Format 0 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)1.000000:0.642905)1.000000:0.567737);
# Format 1 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)E:0.642905)C:0.567737);
# Format 2 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)1.000000:0.642905)1.000000:0.567737);
# Format 3 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)E:0.642905)C:0.567737);
# Format 4 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)));
# Format 5 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729):0.642905):0.567737);
# Format 6 = (A:0.350596,(B:0.728431,(D:0.609498,G:0.125729)E)C);
# Format 7 = (A,(B,(D,G)E)C);
# Format 8 = (A,(B,(D,G)));
# Format 9 = (,(,(,)));

NW_FORMAT = {
    0: [['name', str, True], ["dist", float, True], ['support', float, True], ["dist", float, True]],
    # Flexible with support
    1: [['name', str, True], ["dist", float, True], ['name', str, True], ["dist", float, True]],
    # Flexible with internal node names
    2: [['name', str, False], ["dist", float, False], ['support', float, False], ["dist", float, False]],
    # Strict with support values
    3: [['name', str, False], ["dist", float, False], ['name', str, False], ["dist", float, False]],
    # Strict with internal node names
    4: [['name', str, False], ["dist", float, False], [None, None, False], [None, None, False]],
    5: [['name', str, False], ["dist", float, False], [None, None, False], ["dist", float, False]],
    6: [['name', str, False], [None, None, False], [None, None, False], ["dist", float, False]],
    7: [['name', str, False], ["dist", float, False], ["name", str, False], [None, None, False]],
    8: [['name', str, False], [None, None, False], ["name", str, False], [None, None, False]],
    9: [['name', str, False], [None, None, False], [None, None, False], [None, None, False]],
    # Only topology with node names
    100: [[None, None, False], [None, None, False], [None, None, False], [None, None, False]]  # Only Topology
}


def format_node(node, node_type, format, dist_formatter=None,
                support_formatter=None, name_formatter=None,
                quoted_names=False):
    if dist_formatter is None: dist_formatter = FLOAT_FORMATTER
    if support_formatter is None: support_formatter = FLOAT_FORMATTER
    if name_formatter is None: name_formatter = NAME_FORMATTER

    if node_type == "leaf":
        container1 = NW_FORMAT[format][0][0]  # name
        container2 = NW_FORMAT[format][1][0]  # dists
        converterFn1 = NW_FORMAT[format][0][1]
        converterFn2 = NW_FORMAT[format][1][1]
        flexible1 = NW_FORMAT[format][0][2]
    else:
        container1 = NW_FORMAT[format][2][0]  # support/name
        container2 = NW_FORMAT[format][3][0]  # dist
        converterFn1 = NW_FORMAT[format][2][1]
        converterFn2 = NW_FORMAT[format][3][1]
        flexible1 = NW_FORMAT[format][2][2]

    if converterFn1 == str:
        try:
            if not quoted_names:
                FIRST_PART = re.sub("[" + _ILEGAL_NEWICK_CHARS + "]", "_", \
                                    str(getattr(node, container1)))
            else:
                FIRST_PART = str(getattr(node, container1))
            if not FIRST_PART and container1 == 'name' and not flexible1:
                FIRST_PART = "NoName"

        except (AttributeError, TypeError):
            FIRST_PART = "?"

        FIRST_PART = name_formatter % FIRST_PART
        if quoted_names:
            # FIRST_PART = '"%s"' %FIRST_PART.decode('string_escape').replace('"', '\\"')
            FIRST_PART = '"%s"' % FIRST_PART

    elif converterFn1 is None:
        FIRST_PART = ""
    else:
        try:
            FIRST_PART = support_formatter % (converterFn2(getattr(node, container1)))
        except (ValueError, TypeError):
            FIRST_PART = "?"

    if converterFn2 == str:
        try:
            SECOND_PART = ":" + re.sub("[" + _ILEGAL_NEWICK_CHARS + "]", "_", \
                                       str(getattr(node, container2)))
        except (ValueError, TypeError):
            SECOND_PART = ":?"
    elif converterFn2 is None:
        SECOND_PART = ""
    else:
        try:
            # SECOND_PART = ":%0.6f" %(converterFn2(getattr(node, container2)))
            SECOND_PART = ":%s" % (dist_formatter % (converterFn2(getattr(node, container2))))
        except (ValueError, TypeError):
            SECOND_PART = ":?"

    return "%s%s" % (FIRST_PART, SECOND_PART)


def print_supported_formats():
    from ..coretype.tree import TreeNode
    t = TreeNode()
    t.populate(4, "ABCDEFGHI")
    print(t)
    for f in NW_FORMAT:
        print("Format", f, "=", write_newick(t, features=None, format=f))


class NewickError(Exception):
    """Exception class designed for NewickIO errors."""

    def __init__(self, value):
        if value is None:
            value = ''
        value += "\nYou may want to check other newick loading flags like 'format' or 'quoted_node_names'."
        Exception.__init__(self, value)


def read_newick(newick, root_node=None, format=0, quoted_names=False):
    """ Reads a newick tree from either a string or a file, and returns
    an ETE tree structure.

    A previously existent node object can be passed as the root of the
    tree, which means that all its new children will belong to the same
    class as the root(This allows to work with custom TreeNode
    objects).

    You can also take advantage from this behaviour to concatenate
    several tree structures.
    """

    if root_node is None:
        from ..coretype.tree import TreeNode
        root_node = TreeNode()

    if isinstance(newick, str):
        if os.path.exists(newick):
            if newick.endswith('.gz'):
                import gzip
                nw = gzip.open(newick).read()
            else:
                nw = open(newick, 'rU').read()
        else:
            nw = newick

        matcher = compile_matchers(formatcode=format)
        nw = nw.strip()
        if not nw.startswith('(') and nw.endswith(';'):
            # return _read_node_data(nw[:-1], root_node, "single", matcher, format)
            return _read_newick_from_string(nw, root_node, matcher, format, quoted_names)
        elif not nw.startswith('(') or not nw.endswith(';'):
            raise NewickError('Unexisting tree file or Malformed newick tree structure.')
        else:
            return _read_newick_from_string(nw, root_node, matcher, format, quoted_names)

    else:
        raise NewickError("'newick' argument must be either a filename or a newick string.")


def _read_newick_from_string(nw, root_node, matcher, formatcode, quoted_names):
    """ Reads a newick string in the New Hampshire format. """

    if quoted_names:
        # Quoted text is mapped to references
        quoted_map = {}
        unquoted_nw = ''
        counter = 0
        for token in re.split(_QUOTED_TEXT_RE, nw):
            counter += 1
            if counter % 2 == 1:  # normal newick tree structure data
                unquoted_nw += token
            else:  # quoted text, add to dictionary and replace with reference
                quoted_ref_id = _QUOTED_TEXT_PREFIX + str(int(counter / 2))
                unquoted_nw += quoted_ref_id
                quoted_map[quoted_ref_id] = token[1:-1]  # without the quotes
        nw = unquoted_nw

    if not nw.startswith('(') and nw.endswith(';'):
        _read_node_data(nw[:-1], root_node, "single", matcher, format)
        if quoted_names:
            if root_node.name.startswith(_QUOTED_TEXT_PREFIX):
                root_node.name = quoted_map[root_node.name]
        return root_node

    if nw.count('(') != nw.count(')'):
        raise NewickError('Parentheses do not match. Broken tree structure?')

    # white spaces and separators are removed
    nw = re.sub("[\n\r\t]+", "", nw)

    current_parent = None
    # Each chunk represents the content of a parent node, and it could contain
    # leaves and closing parentheses.
    # We may find:
    # leaf, ..., leaf,
    # leaf, ..., leaf))),
    # leaf)), leaf, leaf))
    # leaf))
    # ) only if formatcode == 100

    for chunk in nw.split("(")[1:]:
        # If no node has been created so far, this is the root, so use the node.
        current_parent = root_node if current_parent is None else current_parent.add_child()

        subchunks = [ch.strip() for ch in chunk.split(",")]
        # We should expect that the chunk finished with a comma (if next chunk
        # is an internal sister node) or a subchunk containing closing parenthesis until the end of the tree.
        # [leaf, leaf, '']
        # [leaf, leaf, ')))', leaf, leaf, '']
        # [leaf, leaf, ')))', leaf, leaf, '']
        # [leaf, leaf, ')))', leaf), leaf, 'leaf);']
        if subchunks[-1] != '' and not subchunks[-1].endswith(';'):
            raise NewickError('Broken newick structure at: %s' % chunk)

        # lets process the subchunks. Every closing parenthesis will close a
        # node and go up one level.
        for i, leaf in enumerate(subchunks):
            if leaf.strip() == '' and i == len(subchunks) - 1:
                continue  # "blah blah ,( blah blah"
            closing_nodes = leaf.split(")")

            # first part after splitting by ) always contain leaf info
            _read_node_data(closing_nodes[0], current_parent, "leaf", matcher, formatcode)

            # next contain closing nodes and data about the internal nodes.
            if len(closing_nodes) > 1:
                for closing_internal in closing_nodes[1:]:
                    closing_internal = closing_internal.rstrip(";")
                    # read internal node data and go up one level
                    _read_node_data(closing_internal, current_parent, "internal", matcher, formatcode)
                    current_parent = current_parent.up

    # references in node names are replaced with quoted text before returning
    if quoted_names:
        for node in root_node.traverse():
            if node.name.startswith(_QUOTED_TEXT_PREFIX):
                node.name = quoted_map[node.name]

    return root_node


def _parse_extra_features(node, NHX_string):
    """ Reads node's extra data form its NHX string. NHX uses this
    format:  [&&NHX:prop1=value1:prop2=value2] """
    NHX_string = NHX_string.replace("[&&NHX:", "")
    NHX_string = NHX_string.replace("]", "")
    for field in NHX_string.split(":"):
        try:
            pname, pvalue = field.split("=")
        except ValueError as e:
            raise NewickError('Invalid NHX format %s' % field)
        node.add_feature(pname, pvalue)


def compile_matchers(formatcode):
    matchers = {}
    for node_type in ["leaf", "single", "internal"]:
        if node_type == "leaf" or node_type == "single":
            container1 = NW_FORMAT[formatcode][0][0]
            container2 = NW_FORMAT[formatcode][1][0]
            converterFn1 = NW_FORMAT[formatcode][0][1]
            converterFn2 = NW_FORMAT[formatcode][1][1]
            flexible1 = NW_FORMAT[formatcode][0][2]
            flexible2 = NW_FORMAT[formatcode][1][2]
        else:
            container1 = NW_FORMAT[formatcode][2][0]
            container2 = NW_FORMAT[formatcode][3][0]
            converterFn1 = NW_FORMAT[formatcode][2][1]
            converterFn2 = NW_FORMAT[formatcode][3][1]
            flexible1 = NW_FORMAT[formatcode][2][2]
            flexible2 = NW_FORMAT[formatcode][3][2]

        if converterFn1 == str:
            FIRST_MATCH = "(" + _NAME_RE + ")"
        elif converterFn1 == float:
            FIRST_MATCH = "(" + _FLOAT_RE + ")"
        elif converterFn1 is None:
            FIRST_MATCH = '()'

        if converterFn2 == str:
            SECOND_MATCH = "(:" + _NAME_RE + ")"
        elif converterFn2 == float:
            SECOND_MATCH = "(:" + _FLOAT_RE + ")"
        elif converterFn2 is None:
            SECOND_MATCH = '()'

        if flexible1 and node_type != 'leaf':
            FIRST_MATCH += "?"
        if flexible2:
            SECOND_MATCH += "?"

        matcher_str = '^\s*%s\s*%s\s*(%s)?\s*$' % (FIRST_MATCH, SECOND_MATCH, _NHX_RE)
        compiled_matcher = re.compile(matcher_str)
        matchers[node_type] = [container1, container2, converterFn1, converterFn2, compiled_matcher]

    return matchers


def _read_node_data(subnw, current_node, node_type, matcher, formatcode):
    """ Reads a leaf node from a subpart of the original newick
    tree """

    if node_type == "leaf" or node_type == "single":
        if node_type == "leaf":
            node = current_node.add_child()
        else:
            node = current_node
    else:
        node = current_node

    subnw = subnw.strip()

    if not subnw and node_type == 'leaf' and formatcode != 100:
        raise NewickError('Empty leaf node found')
    elif not subnw:
        return

    container1, container2, converterFn1, converterFn2, compiled_matcher = matcher[node_type]
    data = re.match(compiled_matcher, subnw)
    if data:
        data = data.groups()
        # This prevents ignoring errors even in flexible nodes:
        if subnw and data[0] is None and data[1] is None and data[2] is None:
            raise NewickError("Unexpected newick format '%s'" % subnw)

        if data[0] is not None and data[0] != '':
            node.add_feature(container1, converterFn1(data[0].strip()))

        if data[1] is not None and data[1] != '':
            node.add_feature(container2, converterFn2(data[1][1:].strip()))

        if data[2] is not None \
                and data[2].startswith("[&&NHX"):
            _parse_extra_features(node, data[2])
    else:
        raise NewickError("Unexpected newick format '%s' " % subnw[0:50])
    return


def write_newick(rootnode, features=None, format=1, format_root_node=True,
                 is_leaf_fn=None, dist_formatter=None, support_formatter=None,
                 name_formatter=None, quoted_names=False):
    """ Iteratively export a tree structure and returns its NHX
    representation. """
    newick = []
    leaf = is_leaf_fn if is_leaf_fn else lambda n: not bool(n.children)
    for postorder, node in rootnode.iter_prepostorder(is_leaf_fn=is_leaf_fn):
        if postorder:
            newick.append(")")
            if node.up is not None or format_root_node:
                newick.append(format_node(node, "internal", format,
                                          dist_formatter=dist_formatter,
                                          support_formatter=support_formatter,
                                          name_formatter=name_formatter,
                                          quoted_names=quoted_names))
                newick.append(_get_features_string(node, features))
        else:
            if node is not rootnode and node != node.up.children[0]:
                newick.append(",")

            if leaf(node):
                newick.append(format_node(node, "leaf", format,
                                          dist_formatter=dist_formatter,
                                          support_formatter=support_formatter,
                                          name_formatter=name_formatter,
                                          quoted_names=quoted_names))
                newick.append(_get_features_string(node, features))
            else:
                newick.append("(")

    newick.append(";")
    return ''.join(newick)


def _get_features_string(self, features=None):
    """ Generates the extended newick string NHX with extra data about
    a node. """
    string = ""
    if features is None:
        features = []
    elif features == []:
        features = sorted(self.features)

    for pr in features:
        if hasattr(self, pr):
            raw = getattr(self, pr)
            if type(raw) in ITERABLE_TYPES:
                raw = '|'.join(map(str, raw))
            elif type(raw) == dict:
                raw = '|'.join(map(lambda x, y: "%s-%s" % (x, y), raw.items()))
            elif type(raw) == str:
                pass
            else:
                raw = str(raw)

            value = re.sub("[" + _ILEGAL_NEWICK_CHARS + "]", "_", \
                           raw)
            if string != "":
                string += ":"
            string += "%s=%s" % (pr, str(value))
    if string != "":
        string = "[&&NHX:" + string + "]"

    return string

####------------------------------------------------------######
#### This is from ete3 coretype/tree.py


DEFAULT_COMPACT = False
DEFAULT_SHOWINTERNAL = False
DEFAULT_DIST = 1.0
DEFAULT_SUPPORT = 1.0
DEFAULT_NAME = ""





class TreeError(Exception):
    """
    A problem occurred during a TreeNode operation
    """

    def __init__(self, value=''):
        self.value = value

    def __str__(self):
        return repr(self.value)


class TreeNode(object):
    """
    TreeNode (Tree) class is used to store a tree structure. A tree
    consists of a collection of TreeNode instances connected in a
    hierarchical way. Trees can be loaded from the New Hampshire Newick
    format (newick).

    :argument newick: Path to the file containing the tree or, alternatively,
       the text string containing the same information.

    :argument 0 format: subnewick format

      .. table::

          ======  ==============================================
          FORMAT  DESCRIPTION
          ======  ==============================================
          0        flexible with support values
          1        flexible with internal node names
          2        all branches + leaf names + internal supports
          3        all branches + all names
          4        leaf branches + leaf names
          5        internal and leaf branches + leaf names
          6        internal branches + leaf names
          7        leaf branches + all names
          8        all names
          9        leaf names
          100      topology only
          ======  ==============================================

    :returns: a tree node object which represents the base of the tree.

    **Examples:**

    ::

        t1 = Tree() # creates an empty tree
        t2 = Tree('(A:1,(B:1,(C:1,D:1):0.5):0.5);')
        t3 = Tree('/home/user/myNewickFile.txt')
    """

    def _get_dist(self):
        return self._dist

    def _set_dist(self, value):
        try:
            self._dist = float(value)
        except ValueError:
            raise TreeError('node dist must be a float number')

    def _get_support(self):
        return self._support

    def _set_support(self, value):
        try:
            self._support = float(value)
        except ValueError:
            raise TreeError('node support must be a float number')

    def _get_up(self):
        return self._up

    def _set_up(self, value):
        if type(value) == type(self) or value is None:
            self._up = value
        else:
            raise TreeError("bad node_up type")

    def _get_children(self):
        return self._children

    def _set_children(self, value):
        if type(value) == list and \
                len(set([type(n) == type(self) for n in value])) < 2:
            self._children = value
        else:
            raise TreeError("Incorrect children type")

    def _get_style(self):
        if self._img_style is None:
            self._set_style(None)

        return self._img_style

    def _set_style(self, value):
        self.set_style(value)

    #: Branch length distance to parent node. Default = 0.0
    img_style = property(fget=_get_style, fset=_set_style)

    #: Branch length distance to parent node. Default = 0.0
    dist = property(fget=_get_dist, fset=_set_dist)
    #: Branch support for current node
    support = property(fget=_get_support, fset=_set_support)
    #: Pointer to parent node
    up = property(fget=_get_up, fset=_set_up)
    #: A list of children nodes
    children = property(fget=_get_children, fset=_set_children)

    def _set_face_areas(self, value):
        if isinstance(value, _FaceAreas):
            self._faces = value
        else:
            raise ValueError("[%s] is not a valid FaceAreas instance" % type(value))

    def _get_face_areas(self):
        if not hasattr(self, "_faces"):
            self._faces = _FaceAreas()
        return self._faces

    faces = property(fget=_get_face_areas, \
                     fset=_set_face_areas)

    def __init__(self, newick=None, format=0, dist=None, support=None,
                 name=None, quoted_node_names=False):
        self._children = []
        self._up = None
        self._dist = DEFAULT_DIST
        self._support = DEFAULT_SUPPORT
        self._img_style = None
        self.features = set([])
        # Add basic features
        self.features.update(["dist", "support", "name"])
        if dist is not None:
            self.dist = dist
        if support is not None:
            self.support = support

        self.name = name if name is not None else DEFAULT_NAME

        # Initialize tree
        if newick is not None:
            self._dist = 0.0
            read_newick(newick, root_node=self, format=format,
                        quoted_names=quoted_node_names)

    def __nonzero__(self):
        return True

    def __bool__(self):
        """
        Python3's equivalent of __nonzero__
        If this is not defined bool(class_instance) will call
        __len__ in python3
        """
        return True

    def __repr__(self):
        return "Tree node '%s' (%s)" % (self.name, hex(self.__hash__()))

    def __and__(self, value):
        """ This allows to execute tree&'A' to obtain the descendant node
        whose name is A"""
        value = str(value)
        try:
            first_match = next(self.iter_search_nodes(name=value))
            return first_match
        except StopIteration:
            raise TreeError("Node not found")

    def __add__(self, value):
        """ This allows to sum two trees."""
        # Should a make the sum with two copies of the original trees?
        if type(value) == self.__class__:
            new_root = self.__class__()
            new_root.add_child(self)
            new_root.add_child(value)
            return new_root
        else:
            raise TreeError("Invalid node type")

    def __str__(self):
        """ Print tree in newick format. """
        return self.get_ascii(compact=DEFAULT_COMPACT, \
                              show_internal=DEFAULT_SHOWINTERNAL)

    def __contains__(self, item):
        """ Check if item belongs to this node. The 'item' argument must
        be a node instance or its associated name."""
        if isinstance(item, self.__class__):
            return item in set(self.get_descendants())
        elif type(item) == str:
            return item in set([n.name for n in self.traverse()])

    def __len__(self):
        """Node len returns number of children."""
        return len(self.get_leaves())

    def __iter__(self):
        """ Iterator over leaf nodes"""
        return self.iter_leaves()

    def add_feature(self, pr_name, pr_value):
        """
        Add or update a node's feature.
        """
        setattr(self, pr_name, pr_value)
        self.features.add(pr_name)

    def add_features(self, **features):
        """
        Add or update several features. """
        for fname, fvalue in features.items():
            setattr(self, fname, fvalue)
            self.features.add(fname)

    def del_feature(self, pr_name):
        """
        Permanently deletes a node's feature.
        """
        if hasattr(self, pr_name):
            delattr(self, pr_name)
            self.features.remove(pr_name)

    # Topology management
    def add_child(self, child=None, name=None, dist=None, support=None):
        """
        Adds a new child to this node. If child node is not suplied
        as an argument, a new node instance will be created.

        :argument None child: the node instance to be added as a child.
        :argument None name: the name that will be given to the child.
        :argument None dist: the distance from the node to the child.
        :argument None support: the support value of child partition.

        :returns: The child node instance

        """
        if child is None:
            child = self.__class__()

        if name is not None:
            child.name = name
        if dist is not None:
            child.dist = dist
        if support is not None:
            child.support = support

        self.children.append(child)
        child.up = self
        return child

    def remove_child(self, child):
        """
        Removes a child from this node (parent and child
        nodes still exit but are no longer connected).
        """
        try:
            self.children.remove(child)
        except ValueError as e:
            raise TreeError("child not found")
        else:
            child.up = None
            return child

    def add_sister(self, sister=None, name=None, dist=None):
        """
        Adds a sister to this node. If sister node is not supplied
        as an argument, a new TreeNode instance will be created and
        returned.
        """
        if self.up is None:
            raise TreeError("A parent node is required to add a sister")
        else:
            return self.up.add_child(child=sister, name=name, dist=dist)

    def remove_sister(self, sister=None):
        """
        Removes a sister node. It has the same effect as
        **`TreeNode.up.remove_child(sister)`**

        If a sister node is not supplied, the first sister will be deleted
        and returned.

        :argument sister: A node instance

        :return: The node removed
        """
        sisters = self.get_sisters()
        if len(sisters) > 0:
            if sister is None:
                sister = sisters.pop(0)
            return self.up.remove_child(sister)

    def delete(self, prevent_nondicotomic=True, preserve_branch_length=False):
        """
        Deletes node from the tree structure. Notice that this method
        makes 'disappear' the node from the tree structure. This means
        that children from the deleted node are transferred to the
        next available parent.

        :param True prevent_nondicotomic: When True (default), delete
            function will be execute recursively to prevent
            single-child nodes.

        :param False preserve_branch_length: If True, branch lengths
            of the deleted nodes are transferred (summed up) to its
            parent's branch, thus keeping original distances among
            nodes.

        **Example:**

        ::

                / C
          root-|
               |        / B
                \--- H |
                        \ A

          > H.delete() will produce this structure:

                / C
               |
          root-|--B
               |
                \ A

        """
        parent = self.up
        if parent:
            if preserve_branch_length:
                if len(self.children) == 1:
                    self.children[0].dist += self.dist
                elif len(self.children) > 1:
                    parent.dist += self.dist

            for ch in self.children:
                parent.add_child(ch)

            parent.remove_child(self)

        # Avoids parents with only one child
        if prevent_nondicotomic and parent and \
                len(parent.children) < 2:
            parent.delete(prevent_nondicotomic=False,
                          preserve_branch_length=preserve_branch_length)

    def detach(self):
        """
        Detachs this node (and all its descendants) from its parent
        and returns the referent to itself.

        Detached node conserves all its structure of descendants, and can
        be attached to another node through the 'add_child' function. This
        mechanism can be seen as a cut and paste.
        """

        if self.up:
            self.up.children.remove(self)
            self.up = None
        return self

    def prune(self, nodes, preserve_branch_length=False):
        """Prunes the topology of a node to conserve only the selected list of leaf
        internal nodes. The minimum number of nodes that conserve the
        topological relationships among the requested nodes will be
        retained. Root node is always conserved.

        :var nodes: a list of node names or node objects that should be retained

        :param False preserve_branch_length: If True, branch lengths
          of the deleted nodes are transferred (summed up) to its
          parent's branch, thus keeping original distances among
          nodes.

        **Examples:**

        ::

          t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', format=1)
          t1.prune(['A', 'B'])


          #                /-A
          #          /D /C|
          #       /F|      \-B
          #      |  |
          #    /H|   \-E
          #   |  |                        /-A
          #-root  \-G                 -root
          #   |                           \-B
          #   |   /-I
          #    \K|
          #       \-J



          t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', format=1)
          t1.prune(['A', 'B', 'C'])

          #                /-A
          #          /D /C|
          #       /F|      \-B
          #      |  |
          #    /H|   \-E
          #   |  |                              /-A
          #-root  \-G                  -root- C|
          #   |                                 \-B
          #   |   /-I
          #    \K|
          #       \-J



          t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', format=1)
          t1.prune(['A', 'B', 'I'])


          #                /-A
          #          /D /C|
          #       /F|      \-B
          #      |  |
          #    /H|   \-E                    /-I
          #   |  |                      -root
          #-root  \-G                      |   /-A
          #   |                             \C|
          #   |   /-I                          \-B
          #    \K|
          #       \-J

          t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', format=1)
          t1.prune(['A', 'B', 'F', 'H'])

          #                /-A
          #          /D /C|
          #       /F|      \-B
          #      |  |
          #    /H|   \-E
          #   |  |                              /-A
          #-root  \-G                -root-H /F|
          #   |                                 \-B
          #   |   /-I
          #    \K|
          #       \-J

        """

        def cmp_nodes(x, y):
            # if several nodes are in the same path of two kept nodes,
            # only one should be maintained. This prioritize internal
            # nodes that are already in the to_keep list and then
            # deeper nodes (closer to the leaves).
            if n2depth[x] > n2depth[y]:
                return -1
            elif n2depth[x] < n2depth[y]:
                return 1
            else:
                return 0

        to_keep = set(_translate_nodes(self, *nodes))
        start, node2path = self.get_common_ancestor(to_keep, get_path=True)
        to_keep.add(self)

        # Calculate which kept nodes are visiting the same nodes in
        # their path to the common ancestor.
        n2count = {}
        n2depth = {}
        for seed, path in node2path.items():
            for visited_node in path:
                if visited_node not in n2depth:
                    depth = visited_node.get_distance(start, topology_only=True)
                    n2depth[visited_node] = depth
                if visited_node is not seed:
                    n2count.setdefault(visited_node, set()).add(seed)

        # if several internal nodes are in the path of exactly the same kept
        # nodes, only one (the deepest) should be maintain.
        visitors2nodes = {}
        for node, visitors in n2count.items():
            # keep nodes connection at least two other nodes
            if len(visitors) > 1:
                visitor_key = frozenset(visitors)
                visitors2nodes.setdefault(visitor_key, set()).add(node)

        for visitors, nodes in visitors2nodes.items():
            if not (to_keep & nodes):
                sorted_nodes = sorted(nodes, key=cmp_to_key(cmp_nodes))
                to_keep.add(sorted_nodes[0])

        for n in self.get_descendants('postorder'):
            if n not in to_keep:
                if preserve_branch_length:
                    if len(n.children) == 1:
                        n.children[0].dist += n.dist
                    elif len(n.children) > 1 and n.up:
                        n.up.dist += n.dist

                n.delete(prevent_nondicotomic=False)

    def swap_children(self):
        """
        Swaps current children order.
        """
        if len(self.children) > 1:
            self.children.reverse()

    # #####################
    # Tree traversing
    # #####################

    def get_children(self):
        """
        Returns an independent list of node's children.
        """
        return [ch for ch in self.children]

    def get_sisters(self):
        """
        Returns an independent list of sister nodes.
        """
        if self.up is not None:
            return [ch for ch in self.up.children if ch != self]
        else:
            return []

    def iter_leaves(self, is_leaf_fn=None):
        """
        Returns an iterator over the leaves under this node.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        for n in self.traverse(strategy="preorder", is_leaf_fn=is_leaf_fn):
            if not is_leaf_fn:
                if n.is_leaf():
                    yield n
            else:
                if is_leaf_fn(n):
                    yield n

    def get_leaves(self, is_leaf_fn=None):
        """
        Returns the list of terminal nodes (leaves) under this node.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        return [n for n in self.iter_leaves(is_leaf_fn=is_leaf_fn)]

    def iter_leaf_names(self, is_leaf_fn=None):
        """
        Returns an iterator over the leaf names under this node.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        for n in self.iter_leaves(is_leaf_fn=is_leaf_fn):
            yield n.name

    def get_leaf_names(self, is_leaf_fn=None):
        """
        Returns the list of terminal node names under the current
        node.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        return [name for name in self.iter_leaf_names(is_leaf_fn=is_leaf_fn)]

    def iter_descendants(self, strategy="levelorder", is_leaf_fn=None):
        """
        Returns an iterator over all descendant nodes.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        for n in self.traverse(strategy=strategy, is_leaf_fn=is_leaf_fn):
            if n is not self:
                yield n

    def get_descendants(self, strategy="levelorder", is_leaf_fn=None):
        """
        Returns a list of all (leaves and internal) descendant nodes.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        return [n for n in self.iter_descendants(strategy=strategy, \
                                                 is_leaf_fn=is_leaf_fn)]

    def traverse(self, strategy="levelorder", is_leaf_fn=None):
        """
        Returns an iterator to traverse the tree structure under this
        node.

        :argument "levelorder" strategy: set the way in which tree
           will be traversed. Possible values are: "preorder" (first
           parent and then children) 'postorder' (first children and
           the parent) and "levelorder" (nodes are visited in order
           from root to leaves)

        :argument None is_leaf_fn: If supplied, ``is_leaf_fn``
           function will be used to interrogate nodes about if they
           are terminal or internal. ``is_leaf_fn`` function should
           receive a node instance as first argument and return True
           or False. Use this argument to traverse a tree by
           dynamically collapsing internal nodes matching
           ``is_leaf_fn``.
        """
        if strategy == "preorder":
            return self._iter_descendants_preorder(is_leaf_fn=is_leaf_fn)
        elif strategy == "levelorder":
            return self._iter_descendants_levelorder(is_leaf_fn=is_leaf_fn)
        elif strategy == "postorder":
            return self._iter_descendants_postorder(is_leaf_fn=is_leaf_fn)

    def iter_prepostorder(self, is_leaf_fn=None):
        """
        Iterate over all nodes in a tree yielding every node in both
        pre and post order. Each iteration returns a postorder flag
        (True if node is being visited in postorder) and a node
        instance.
        """
        to_visit = [self]
        if is_leaf_fn is not None:
            _leaf = is_leaf_fn
        else:
            _leaf = self.__class__.is_leaf

        while to_visit:
            node = to_visit.pop(-1)
            try:
                node = node[1]
            except TypeError:
                # PREORDER ACTIONS
                yield (False, node)
                if not _leaf(node):
                    # ADD CHILDREN
                    to_visit.extend(reversed(node.children + [[1, node]]))
            else:
                # POSTORDER ACTIONS
                yield (True, node)

    def _iter_descendants_postorder(self, is_leaf_fn=None):
        to_visit = [self]
        if is_leaf_fn is not None:
            _leaf = is_leaf_fn
        else:
            _leaf = self.__class__.is_leaf

        while to_visit:
            node = to_visit.pop(-1)
            try:
                node = node[1]
            except TypeError:
                # PREORDER ACTIONS
                if not _leaf(node):
                    # ADD CHILDREN
                    to_visit.extend(reversed(node.children + [[1, node]]))
                else:
                    yield node
            else:
                # POSTORDER ACTIONS
                yield node

    def _iter_descendants_levelorder(self, is_leaf_fn=None):
        """
        Iterate over all desdecendant nodes.
        """
        tovisit = deque([self])
        while len(tovisit) > 0:
            node = tovisit.popleft()
            yield node
            if not is_leaf_fn or not is_leaf_fn(node):
                tovisit.extend(node.children)

    def _iter_descendants_preorder(self, is_leaf_fn=None):
        """
        Iterator over all descendant nodes.
        """
        to_visit = deque()
        node = self
        while node is not None:
            yield node
            if not is_leaf_fn or not is_leaf_fn(node):
                to_visit.extendleft(reversed(node.children))
            try:
                node = to_visit.popleft()
            except:
                node = None

    def iter_ancestors(self):
        '''versionadded: 2.2

        Iterates over the list of all ancestor nodes from current node
        to the current tree root.

        '''
        node = self
        while node.up is not None:
            yield node.up
            node = node.up

    def get_ancestors(self):
        '''versionadded: 2.2

        Returns the list of all ancestor nodes from current node to
        the current tree root.

        '''
        return [n for n in self.iter_ancestors()]

    def describe(self):
        """
        Prints general information about this node and its
        connections.
        """
        if len(self.get_tree_root().children) == 2:
            rooting = "Yes"
        elif len(self.get_tree_root().children) > 2:
            rooting = "No"
        else:
            rooting = "No children"
        max_node, max_dist = self.get_farthest_leaf()
        cached_content = self.get_cached_content()
        print("Number of leaf nodes:\t%d" % len(cached_content[self]))
        print("Total number of nodes:\t%d" % len(cached_content))
        print("Rooted:\t%s" % rooting)
        print("Most distant node:\t%s" % max_node.name)
        print("Max. distance:\t%f" % max_dist)

    def write(self, features=None, outfile=None, format=0, is_leaf_fn=None,
              format_root_node=False, dist_formatter=None, support_formatter=None,
              name_formatter=None, quoted_node_names=False):
        """
        Returns the newick representation of current node. Several
        arguments control the way in which extra data is shown for
        every node:

        :argument features: a list of feature names to be exported
          using the Extended Newick Format (i.e. features=["name",
          "dist"]). Use an empty list to export all available features
          in each node (features=[])

        :argument outfile: writes the output to a given file

        :argument format: defines the newick standard used to encode the
          tree. See tutorial for details.

        :argument False format_root_node: If True, it allows features
          and branch information from root node to be exported as a
          part of the newick text string. For newick compatibility
          reasons, this is False by default.

        :argument is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.

        **Example:**

        ::

             t.get_newick(features=["species","name"], format=1)

        """

        nw = write_newick(self, features=features, format=format,
                          is_leaf_fn=is_leaf_fn,
                          format_root_node=format_root_node,
                          dist_formatter=dist_formatter,
                          support_formatter=support_formatter,
                          name_formatter=name_formatter,
                          quoted_names=quoted_node_names)

        if outfile is not None:
            with open(outfile, "w") as OUT:
                OUT.write(nw)
        else:
            return nw

    def get_tree_root(self):
        """
        Returns the absolute root node of current tree structure.
        """
        root = self
        while root.up is not None:
            root = root.up
        return root

    def get_common_ancestor(self, *target_nodes, **kargs):
        """
        Returns the first common ancestor between this node and a given
        list of 'target_nodes'.

        **Examples:**

        ::

          t = tree.Tree("(((A:0.1, B:0.01):0.001, C:0.0001):1.0[&&NHX:name=common], (D:0.00001):0.000001):2.0[&&NHX:name=root];")
          A = t.get_descendants_by_name("A")[0]
          C = t.get_descendants_by_name("C")[0]
          common =  A.get_common_ancestor(C)
          print common.name

        """

        get_path = kargs.get("get_path", False)

        if len(target_nodes) == 1 and type(target_nodes[0]) \
                in set([set, tuple, list, frozenset]):
            target_nodes = target_nodes[0]

        # Convert node names into node instances
        target_nodes = _translate_nodes(self, *target_nodes)

        # If only one node is provided, use self as the second target
        if type(target_nodes) != list:
            target_nodes = [target_nodes, self]

        n2path = {}
        reference = []
        ref_node = None
        for n in target_nodes:
            current = n
            while current:
                n2path.setdefault(n, set()).add(current)
                if not ref_node:
                    reference.append(current)
                current = current.up
            if not ref_node:
                ref_node = n

        common = None
        for n in reference:
            broken = False
            for node, path in n2path.items():
                if node is not ref_node and n not in path:
                    broken = True
                    break

            if not broken:
                common = n
                break
        if not common:
            raise TreeError("Nodes are not connected!")

        if get_path:
            return common, n2path
        else:
            return common

    def iter_search_nodes(self, **conditions):
        """
        Search nodes in an iterative way. Matches are yielded as they
        are being found. This avoids needing to scan the full tree
        topology before returning the first matches. Useful when
        dealing with huge trees.
        """

        for n in self.traverse():
            conditions_passed = 0
            for key, value in conditions.items():
                if hasattr(n, key) and getattr(n, key) == value:
                    conditions_passed += 1
            if conditions_passed == len(conditions):
                yield n

    def search_nodes(self, **conditions):
        """
        Returns the list of nodes matching a given set of conditions.

        **Example:**

        ::

          tree.search_nodes(dist=0.0, name="human")

        """
        matching_nodes = []
        for n in self.iter_search_nodes(**conditions):
            matching_nodes.append(n)
        return matching_nodes

    def get_leaves_by_name(self, name):
        """
        Returns a list of leaf nodes matching a given name.
        """
        return self.search_nodes(name=name, children=[])

    def is_leaf(self):
        """
        Return True if current node is a leaf.
        """
        return len(self.children) == 0

    def is_root(self):
        """
        Returns True if current node has no parent
        """
        if self.up is None:
            return True
        else:
            return False

    # ###########################
    # Distance related functions
    # ###########################
    def get_distance(self, target, target2=None, topology_only=False):
        """
        Returns the distance between two nodes. If only one target is
        specified, it returns the distance between the target and the
        current node.

        :argument target: a node within the same tree structure.

        :argument target2: a node within the same tree structure. If
          not specified, current node is used as target2.

        :argument False topology_only: If set to True, distance will
          refer to the number of nodes between target and target2.

        :returns: branch length distance between target and
          target2. If topology_only flag is True, returns the number
          of nodes between target and target2.

        """

        if target2 is None:
            target2 = self
            root = self.get_tree_root()
        else:
            # is target node under current node?
            root = self

        target, target2 = _translate_nodes(root, target, target2)
        ancestor = root.get_common_ancestor(target, target2)

        dist = 0.0
        for n in [target2, target]:
            current = n
            while current != ancestor:
                if topology_only:
                    if current != target:
                        dist += 1
                else:
                    dist += current.dist
                current = current.up
        return dist

    def get_farthest_node(self, topology_only=False):
        """
        Returns the node's farthest descendant or ancestor node, and the
        distance to it.

        :argument False topology_only: If set to True, distance
          between nodes will be referred to the number of nodes
          between them. In other words, topological distance will be
          used instead of branch length distances.

        :return: A tuple containing the farthest node referred to the
          current node and the distance to it.

        """
        # Init farthest node to current farthest leaf
        farthest_node, farthest_dist = self.get_farthest_leaf(topology_only=topology_only)

        prev = self
        cdist = 0.0 if topology_only else prev.dist
        current = prev.up
        while current is not None:
            for ch in current.children:
                if ch != prev:
                    if not ch.is_leaf():
                        fnode, fdist = ch.get_farthest_leaf(topology_only=topology_only)
                    else:
                        fnode = ch
                        fdist = 0
                    if topology_only:
                        fdist += 1.0
                    else:
                        fdist += ch.dist
                    if cdist + fdist > farthest_dist:
                        farthest_dist = cdist + fdist
                        farthest_node = fnode
            prev = current
            if topology_only:
                cdist += 1
            else:
                cdist += prev.dist
            current = prev.up
        return farthest_node, farthest_dist

    def _get_farthest_and_closest_leaves(self, topology_only=False, is_leaf_fn=None):
        # if called from a leaf node, no necessary to compute
        if (is_leaf_fn and is_leaf_fn(self)) or self.is_leaf():
            return self, 0.0, self, 0.0

        min_dist = None
        min_node = None
        max_dist = None
        max_node = None
        d = 0.0
        for post, n in self.iter_prepostorder(is_leaf_fn=is_leaf_fn):
            if n is self:
                continue
            if post:
                d -= n.dist if not topology_only else 1.0
            else:
                if (is_leaf_fn and is_leaf_fn(n)) or n.is_leaf():
                    total_d = d + n.dist if not topology_only else d
                    if min_dist is None or total_d < min_dist:
                        min_dist = total_d
                        min_node = n
                    if max_dist is None or total_d > max_dist:
                        max_dist = total_d
                        max_node = n
                else:
                    d += n.dist if not topology_only else 1.0
        return min_node, min_dist, max_node, max_dist

    def get_farthest_leaf(self, topology_only=False, is_leaf_fn=None):
        """
        Returns node's farthest descendant node (which is always a leaf), and the
        distance to it.

        :argument False topology_only: If set to True, distance
          between nodes will be referred to the number of nodes
          between them. In other words, topological distance will be
          used instead of branch length distances.

        :return: A tuple containing the farthest leaf referred to the
          current node and the distance to it.
        """
        min_node, min_dist, max_node, max_dist = self._get_farthest_and_closest_leaves(
            topology_only=topology_only, is_leaf_fn=is_leaf_fn)
        return max_node, max_dist

    def get_closest_leaf(self, topology_only=False, is_leaf_fn=None):
        """Returns node's closest descendant leaf and the distance to
        it.

        :argument False topology_only: If set to True, distance
          between nodes will be referred to the number of nodes
          between them. In other words, topological distance will be
          used instead of branch length distances.

        :return: A tuple containing the closest leaf referred to the
          current node and the distance to it.

        """
        min_node, min_dist, max_node, max_dist = self._get_farthest_and_closest_leaves(
            topology_only=topology_only, is_leaf_fn=is_leaf_fn)

        return min_node, min_dist

    def get_midpoint_outgroup(self):
        """
        Returns the node that divides the current tree into two distance-balanced
        partitions.
        """
        # Gets the farthest node to the current root
        root = self.get_tree_root()
        nA, r2A_dist = root.get_farthest_leaf()
        nB, A2B_dist = nA.get_farthest_node()

        outgroup = nA
        middist = A2B_dist / 2.0
        cdist = 0
        current = nA
        while current is not None:
            cdist += current.dist
            if cdist > (middist):  # Deja de subir cuando se pasa del maximo
                break
            else:
                current = current.up
        return current

    def populate(self, size, names_library=None, reuse_names=False,
                 random_branches=False, branch_range=(0, 1),
                 support_range=(0, 1)):
        """
        Generates a random topology by populating current node.

        :argument None names_library: If provided, names library
          (list, set, dict, etc.) will be used to name nodes.

        :argument False reuse_names: If True, node names will not be
          necessarily unique, which makes the process a bit more
          efficient.

        :argument False random_branches: If True, branch distances and support
          values will be randomized.

        :argument (0,1) branch_range: If random_branches is True, this
          range of values will be used to generate random distances.

        :argument (0,1) support_range: If random_branches is True,
          this range of values will be used to generate random branch
          support values.

        """
        NewNode = self.__class__

        if len(self.children) > 1:
            connector = NewNode()
            for ch in self.get_children():
                ch.detach()
                connector.add_child(child=ch)
            root = NewNode()
            self.add_child(child=connector)
            self.add_child(child=root)
        else:
            root = self

        next_deq = deque([root])
        for i in range(size - 1):
            if random.randint(0, 1):
                p = next_deq.pop()
            else:
                p = next_deq.popleft()

            c1 = p.add_child()
            c2 = p.add_child()
            next_deq.extend([c1, c2])
            if random_branches:
                c1.dist = random.uniform(*branch_range)
                c2.dist = random.uniform(*branch_range)
                c1.support = random.uniform(*branch_range)
                c2.support = random.uniform(*branch_range)
            else:
                c1.dist = 1.0
                c2.dist = 1.0
                c1.support = 1.0
                c2.support = 1.0

        # next contains leaf nodes
        charset = "abcdefghijklmnopqrstuvwxyz"
        if names_library:
            names_library = deque(names_library)
        else:
            avail_names = itertools.combinations_with_replacement(charset, 10)
        for n in next_deq:
            if names_library:
                if reuse_names:
                    tname = random.sample(names_library, 1)[0]
                else:
                    tname = names_library.pop()
            else:
                tname = ''.join(next(avail_names))
            n.name = tname

    def set_outgroup(self, outgroup):
        """
        Sets a descendant node as the outgroup of a tree.  This function
        can be used to root a tree or even an internal node.

        :argument outgroup: a node instance within the same tree
          structure that will be used as a basal node.

        """

        outgroup = _translate_nodes(self, outgroup)

        if self == outgroup:
            raise TreeError("Cannot set myself as outgroup")

        parent_outgroup = outgroup.up

        # Detects (sub)tree root
        n = outgroup
        while n.up is not self:
            n = n.up

        # If outgroup is a child from root, but with more than one
        # sister nodes, creates a new node to group them

        self.children.remove(n)
        if len(self.children) != 1:
            down_branch_connector = self.__class__()
            down_branch_connector.dist = 0.0
            down_branch_connector.support = n.support
            for ch in self.get_children():
                down_branch_connector.children.append(ch)
                ch.up = down_branch_connector
                self.children.remove(ch)
        else:
            down_branch_connector = self.children[0]

        # Connects down branch to myself or to outgroup
        quien_va_ser_padre = parent_outgroup
        if quien_va_ser_padre is not self:
            # Parent-child swapping
            quien_va_ser_hijo = quien_va_ser_padre.up
            quien_fue_padre = None
            buffered_dist = quien_va_ser_padre.dist
            buffered_support = quien_va_ser_padre.support

            while quien_va_ser_hijo is not self:
                quien_va_ser_padre.children.append(quien_va_ser_hijo)
                quien_va_ser_hijo.children.remove(quien_va_ser_padre)

                buffered_dist2 = quien_va_ser_hijo.dist
                buffered_support2 = quien_va_ser_hijo.support
                quien_va_ser_hijo.dist = buffered_dist
                quien_va_ser_hijo.support = buffered_support
                buffered_dist = buffered_dist2
                buffered_support = buffered_support2

                quien_va_ser_padre.up = quien_fue_padre
                quien_fue_padre = quien_va_ser_padre

                quien_va_ser_padre = quien_va_ser_hijo
                quien_va_ser_hijo = quien_va_ser_padre.up

            quien_va_ser_padre.children.append(down_branch_connector)
            down_branch_connector.up = quien_va_ser_padre
            quien_va_ser_padre.up = quien_fue_padre

            down_branch_connector.dist += buffered_dist
            outgroup2 = parent_outgroup
            parent_outgroup.children.remove(outgroup)
            outgroup2.dist = 0

        else:
            outgroup2 = down_branch_connector

        outgroup.up = self
        outgroup2.up = self
        # outgroup is always the first children. Some function my
        # trust on this fact, so do no change this.
        self.children = [outgroup, outgroup2]
        middist = (outgroup2.dist + outgroup.dist) / 2
        outgroup.dist = middist
        outgroup2.dist = middist
        outgroup2.support = outgroup.support

    def unroot(self):
        """
        Unroots current node. This function is expected to be used on
        the absolute tree root node, but it can be also be applied to
        any other internal node. It will convert a split into a
        multifurcation.
        """
        if len(self.children) == 2:
            if not self.children[0].is_leaf():
                self.children[0].delete()
            elif not self.children[1].is_leaf():
                self.children[1].delete()
            else:
                raise TreeError("Cannot unroot a tree with only two leaves")

    #
    #
    # def copy(self, method="cpickle"):
    #     """.. versionadded: 2.1
    #
    #     Returns a copy of the current node.
    #
    #     :var cpickle method: Protocol used to copy the node
    #     structure. The following values are accepted:
    #
    #        - "newick": Tree topology, node names, branch lengths and
    #          branch support values will be copied by as represented in
    #          the newick string (copy by newick string serialisation).
    #
    #        - "newick-extended": Tree topology and all node features
    #          will be copied based on the extended newick format
    #          representation. Only node features will be copied, thus
    #          excluding other node attributes. As this method is also
    #          based on newick serialisation, features will be converted
    #          into text strings when making the copy.
    #
    #        - "cpickle": The whole node structure and its content is
    #          cloned based on cPickle object serialisation (slower, but
    #          recommended for full tree copying)
    #
    #        - "deepcopy": The whole node structure and its content is
    #          copied based on the standard "copy" Python functionality
    #          (this is the slowest method but it allows to copy complex
    #          objects even if attributes point to lambda functions,
    #          etc.)
    #
    #     """
    #     method = method.lower()
    #     if method == "newick":
    #         new_node = self.__class__(self.write(features=["name"], format_root_node=True))
    #     elif method == "newick-extended":
    #         self.write(features=[], format_root_node=True)
    #         new_node = self.__class__(self.write(features=[]))
    #     elif method == "deepcopy":
    #         parent = self.up
    #         self.up = None
    #         new_node = copy.deepcopy(self)
    #         self.up = parent
    #     elif method == "cpickle":
    #         parent = self.up
    #         self.up = None
    #         new_node = six.moves.cPickle.loads(six.moves.cPickle.dumps(self, 2))
    #         self.up = parent
    #     else:
    #         raise TreeError("Invalid copy method")
    #
    #     return new_node

    def _asciiArt(self, char1='-', show_internal=True, compact=False, attributes=None):
        """
        Returns the ASCII representation of the tree.

        Code based on the PyCogent GPL project.
        """
        if not attributes:
            attributes = ["name"]
        node_name = ', '.join(map(str, [getattr(self, v) for v in attributes if hasattr(self, v)]))

        LEN = max(3, len(node_name) if not self.children or show_internal else 3)
        PAD = ' ' * LEN
        PA = ' ' * (LEN - 1)
        if not self.is_leaf():
            mids = []
            result = []
            for c in self.children:
                if len(self.children) == 1:
                    char2 = '/'
                elif c is self.children[0]:
                    char2 = '/'
                elif c is self.children[-1]:
                    char2 = '\\'
                else:
                    char2 = '-'
                (clines, mid) = c._asciiArt(char2, show_internal, compact, attributes)
                mids.append(mid + len(result))
                result.extend(clines)
                if not compact:
                    result.append('')
            if not compact:
                result.pop()
            (lo, hi, end) = (mids[0], mids[-1], len(result))
            prefixes = [PAD] * (lo + 1) + [PA + '|'] * (hi - lo - 1) + [PAD] * (end - hi)
            mid = int((lo + hi) / 2)
            prefixes[mid] = char1 + '-' * (LEN - 2) + prefixes[mid][-1]
            result = [p + l for (p, l) in zip(prefixes, result)]
            if show_internal:
                stem = result[mid]
                result[mid] = stem[0] + node_name + stem[len(node_name) + 1:]
            return (result, mid)
        else:
            return ([char1 + '-' + node_name], 0)

    def get_ascii(self, show_internal=True, compact=False, attributes=None):
        """
        Returns a string containing an ascii drawing of the tree.

        :argument show_internal: includes internal edge names.
        :argument compact: use exactly one line per tip.

        :param attributes: A list of node attributes to shown in the
            ASCII representation.

        """
        (lines, mid) = self._asciiArt(show_internal=show_internal,
                                      compact=compact, attributes=attributes)
        return '\n' + '\n'.join(lines)

    def ladderize(self, direction=0):
        """
        .. versionadded: 2.1

        Sort the branches of a given tree (swapping children nodes)
        according to the size of each partition.

        ::

           t =  Tree("(f,((d, ((a,b),c)),e));")

           print t

           #
           #      /-f
           #     |
           #     |          /-d
           # ----|         |
           #     |     /---|          /-a
           #     |    |    |     /---|
           #     |    |     \---|     \-b
           #      \---|         |
           #          |          \-c
           #          |
           #           \-e

           t.ladderize()
           print t

           #      /-f
           # ----|
           #     |     /-e
           #      \---|
           #          |     /-d
           #           \---|
           #               |     /-c
           #                \---|
           #                    |     /-a
           #                     \---|
           #                          \-b

        """

        if not self.is_leaf():
            n2s = {}
            for n in self.get_children():
                s = n.ladderize(direction=direction)
                n2s[n] = s

            self.children.sort(key=lambda x: n2s[x])
            if direction == 1:
                self.children.reverse()
            size = sum(n2s.values())
        else:
            size = 1

        return size

    def sort_descendants(self, attr="name"):
        """
        .. versionadded: 2.1

        Sort the branches of a given tree by node names. After the
        tree is sorted, nodes are labeled in ascending order. This
        can be used to ensure that nodes in a tree with the same node
        names are always labeled in the same way. Note that if
        duplicated names are present, extra criteria should be added
        to sort nodes.

        Unique id is stored as a node._nid attribute

        """

        node2content = self.get_cached_content(store_attr=attr, container_type=list)

        for n in self.traverse():
            if not n.is_leaf():
                n.children.sort(key=lambda x: str(sorted(node2content[x])))

    def get_cached_content(self, store_attr=None, container_type=set, leaves_only=True, _store=None):
        """
        .. versionadded: 2.2

        Returns a dictionary pointing to the preloaded content of each
        internal node under this tree. Such a dictionary is intended
        to work as a cache for operations that require many traversal
        operations.

        :param None store_attr: Specifies the node attribute that
            should be cached (i.e. name, distance, etc.). When none,
            the whole node instance is cached.

        :param _store: (internal use)

        """

        if _store is None:
            _store = {}

        def get_value(_n):
            if store_attr is None:
                _val = [_n]
            else:
                if not isinstance(store_attr, str):
                    _val = [tuple(getattr(_n, attr, None) for attr in store_attr)]

                else:
                    _val = [getattr(_n, store_attr, None)]

            return _val

        for ch in self.children:
            ch.get_cached_content(store_attr=store_attr,
                                  container_type=container_type,
                                  leaves_only=leaves_only,
                                  _store=_store)

        if self.children:
            if not leaves_only:
                val = container_type(get_value(self))
            else:
                val = container_type()
            for ch in self.children:
                if type(val) == list:
                    val.extend(_store[ch])
                if type(val) == set:
                    val.update(_store[ch])

                if not leaves_only:
                    if type(val) == list:
                        val.extend(get_value(ch))
                    if type(val) == set:
                        val.update(get_value(ch))

            _store[self] = val
        else:
            _store[self] = container_type(get_value(self))

        return _store

    def robinson_foulds(self, t2, attr_t1="name", attr_t2="name",
                        unrooted_trees=False, expand_polytomies=False,
                        polytomy_size_limit=5, skip_large_polytomies=False,
                        correct_by_polytomy_size=False, min_support_t1=0.0,
                        min_support_t2=0.0):
        """
        .. versionadded: 2.2

        Returns the Robinson-Foulds symmetric distance between current
        tree and a different tree instance.

        :param t2: reference tree

        :param name attr_t1: Compare trees using a custom node
                              attribute as a node name.

        :param name attr_t2: Compare trees using a custom node
                              attribute as a node name in target tree.

        :param False attr_t2: If True, consider trees as unrooted.

        :param False expand_polytomies: If True, all polytomies in the reference
           and target tree will be expanded into all possible binary
           trees. Robinson-foulds distance will be calculated between all
           tree combinations and the minimum value will be returned.
           See also, :func:`NodeTree.expand_polytomy`.

        :returns: (rf, rf_max, common_attrs, names, edges_t1, edges_t2,  discarded_edges_t1, discarded_edges_t2)

        """
        ref_t = self
        target_t = t2
        if not unrooted_trees and (len(ref_t.children) > 2 or len(target_t.children) > 2):
            raise TreeError("Unrooted tree found! You may want to activate the unrooted_trees flag.")

        if expand_polytomies and correct_by_polytomy_size:
            raise TreeError("expand_polytomies and correct_by_polytomy_size are mutually exclusive.")

        if expand_polytomies and unrooted_trees:
            raise TreeError("expand_polytomies and unrooted_trees arguments cannot be enabled at the same time")

        attrs_t1 = set([getattr(n, attr_t1) for n in ref_t.iter_leaves() if hasattr(n, attr_t1)])
        attrs_t2 = set([getattr(n, attr_t2) for n in target_t.iter_leaves() if hasattr(n, attr_t2)])
        common_attrs = attrs_t1 & attrs_t2
        # release mem
        attrs_t1, attrs_t2 = None, None

        # Check for duplicated items (is it necessary? can we optimize? what's the impact in performance?')
        size1 = len([True for n in ref_t.iter_leaves() if getattr(n, attr_t1, None) in common_attrs])
        size2 = len([True for n in target_t.iter_leaves() if getattr(n, attr_t2, None) in common_attrs])
        if size1 > len(common_attrs):
            raise TreeError('Duplicated items found in source tree')
        if size2 > len(common_attrs):
            raise TreeError('Duplicated items found in reference tree')

        if expand_polytomies:
            ref_trees = [Tree(nw) for nw in
                         ref_t.expand_polytomies(map_attr=attr_t1,
                                                 polytomy_size_limit=polytomy_size_limit,
                                                 skip_large_polytomies=skip_large_polytomies)]
            target_trees = [Tree(nw) for nw in
                            target_t.expand_polytomies(map_attr=attr_t2,
                                                       polytomy_size_limit=polytomy_size_limit,
                                                       skip_large_polytomies=skip_large_polytomies)]
            attr_t1, attr_t2 = "name", "name"
        else:
            ref_trees = [ref_t]
            target_trees = [target_t]

        polytomy_correction = 0
        if correct_by_polytomy_size:
            corr1 = sum([0] + [len(n.children) - 2 for n in ref_t.traverse() if len(n.children) > 2])
            corr2 = sum([0] + [len(n.children) - 2 for n in target_t.traverse() if len(n.children) > 2])
            if corr1 and corr2:
                raise TreeError("Both trees contain polytomies! Try expand_polytomies=True instead")
            else:
                polytomy_correction = max([corr1, corr2])

        min_comparison = None
        for t1 in ref_trees:
            t1_content = t1.get_cached_content()
            t1_leaves = t1_content[t1]
            if unrooted_trees:
                edges1 = set([
                    tuple(sorted([tuple(sorted([getattr(n, attr_t1) for n in content if
                                                hasattr(n, attr_t1) and getattr(n, attr_t1) in common_attrs])),
                                  tuple(sorted([getattr(n, attr_t1) for n in t1_leaves - content if
                                                hasattr(n, attr_t1) and getattr(n, attr_t1) in common_attrs]))]))
                    for content in t1_content.values()])
                edges1.discard(((), ()))
            else:
                edges1 = set([
                    tuple(sorted([getattr(n, attr_t1) for n in content if
                                  hasattr(n, attr_t1) and getattr(n, attr_t1) in common_attrs]))
                    for content in t1_content.values()])
                edges1.discard(())

            if min_support_t1:
                support_t1 = dict([
                    (tuple(sorted([getattr(n, attr_t1) for n in content if
                                   hasattr(n, attr_t1) and getattr(n, attr_t1) in common_attrs])), branch.support)
                    for branch, content in t1_content.values()])

            for t2 in target_trees:
                t2_content = t2.get_cached_content()
                t2_leaves = t2_content[t2]
                if unrooted_trees:
                    edges2 = set([
                        tuple(sorted([
                            tuple(sorted([getattr(n, attr_t2) for n in content if
                                          hasattr(n, attr_t2) and getattr(n, attr_t2) in common_attrs])),
                            tuple(sorted([getattr(n, attr_t2) for n in t2_leaves - content if
                                          hasattr(n, attr_t2) and getattr(n, attr_t2) in common_attrs]))]))
                        for content in t2_content.values()])
                    edges2.discard(((), ()))
                else:
                    edges2 = set([
                        tuple(sorted([getattr(n, attr_t2) for n in content if
                                      hasattr(n, attr_t2) and getattr(n, attr_t2) in common_attrs]))
                        for content in t2_content.values()])
                    edges2.discard(())

                if min_support_t2:
                    support_t2 = dict([
                        (tuple(sorted(([getattr(n, attr_t2) for n in content if
                                        hasattr(n, attr_t2) and getattr(n, attr_t2) in common_attrs]))), branch.support)
                        for branch, content in t2_content.values()])

                # if a support value is passed as a constraint, discard lowly supported branches from the analysis
                discard_t1, discard_t2 = set(), set()
                if min_support_t1 and unrooted_trees:
                    discard_t1 = set(
                        [p for p in edges1 if support_t1.get(p[0], support_t1.get(p[1], 999999999)) < min_support_t1])
                elif min_support_t1:
                    discard_t1 = set([p for p in edges1 if support_t1[p] < min_support_t1])

                if min_support_t2 and unrooted_trees:
                    discard_t2 = set(
                        [p for p in edges2 if support_t2.get(p[0], support_t2.get(p[1], 999999999)) < min_support_t2])
                elif min_support_t2:
                    discard_t2 = set([p for p in edges2 if support_t2[p] < min_support_t2])

                # rf = len(edges1 ^ edges2) - (len(discard_t1) + len(discard_t2)) - polytomy_correction # poly_corr is 0 if the flag is not enabled
                # rf = len((edges1-discard_t1) ^ (edges2-discard_t2)) - polytomy_correction

                # the two root edges are never counted here, as they are always
                # present in both trees because of the common attr filters
                rf = len(((edges1 ^ edges2) - discard_t2) - discard_t1) - polytomy_correction

                if unrooted_trees:
                    # thought this may work, but it does not, still I don't see why
                    # max_parts = (len(common_attrs)*2) - 6 - len(discard_t1) - len(discard_t2)
                    max_parts = (len([p for p in edges1 - discard_t1 if len(p[0]) > 1 and len(p[1]) > 1]) +
                                 len([p for p in edges2 - discard_t2 if len(p[0]) > 1 and len(p[1]) > 1]))
                else:
                    # thought this may work, but it does not, still I don't see why
                    # max_parts = (len(common_attrs)*2) - 4 - len(discard_t1) - len(discard_t2)

                    # Otherwise we need to count the actual number of valid
                    # partitions in each tree -2 is to avoid counting the root
                    # partition of the two trees (only needed in rooted trees)
                    max_parts = (len([p for p in edges1 - discard_t1 if len(p) > 1]) +
                                 len([p for p in edges2 - discard_t2 if len(p) > 1])) - 2

                    # print max_parts

                if not min_comparison or min_comparison[0] > rf:
                    min_comparison = [rf, max_parts, common_attrs, edges1, edges2, discard_t1, discard_t2]

        return min_comparison

    def compare(self, ref_tree, use_collateral=False, min_support_source=0.0, min_support_ref=0.0,
                has_duplications=False, expand_polytomies=False, unrooted=False,
                max_treeko_splits_to_be_artifact=1000, ref_tree_attr='name', source_tree_attr='name'):

        """compare this tree with another using robinson foulds symmetric difference
        and number of shared edges. Trees of different sizes and with duplicated
        items allowed.

        returns: a Python dictionary with results

        """
        source_tree = self

        def _safe_div(a, b):
            if a != 0:
                return a / float(b)
            else:
                return 0.0

        def _compare(src_tree, ref_tree):
            # calculate partitions and rf distances
            rf, maxrf, common, ref_p, src_p, ref_disc, src_disc = ref_tree.robinson_foulds(src_tree,
                                                                                           expand_polytomies=expand_polytomies,
                                                                                           unrooted_trees=unrooted,
                                                                                           attr_t1=ref_tree_attr,
                                                                                           attr_t2=source_tree_attr,
                                                                                           min_support_t2=min_support_source,
                                                                                           min_support_t1=min_support_ref)

            # if trees share leaves, count their distances
            if len(common) > 0 and src_p and ref_p:
                if unrooted:
                    valid_ref_edges = set([p for p in (ref_p - ref_disc) if len(p[0]) > 1 and len(p[1]) > 0])
                    valid_src_edges = set([p for p in (src_p - src_disc) if len(p[0]) > 1 and len(p[1]) > 0])
                    common_edges = valid_ref_edges & valid_src_edges
                else:

                    valid_ref_edges = set([p for p in (ref_p - ref_disc) if len(p) > 1])
                    valid_src_edges = set([p for p in (src_p - src_disc) if len(p) > 1])
                    common_edges = valid_ref_edges & valid_src_edges

            else:
                valid_ref_edges = set()
                valid_src_edges = set()
                common_edges = set()

                # # % of ref edges found in tree
                # ref_found.append(float(len(p2 & p1)) / reftree_edges)

                # # valid edges in target, discard also leaves
                # p2bis = set([p for p in (p2-d2) if len(p[0])>1 and len(p[1])>1])
                # if p2bis:
                #     incompatible_target_branches = float(len((p2-d2) - p1))
                #     target_found.append(1 - (incompatible_target_branches / (len(p2-d2))))

            return rf, maxrf, len(common), valid_ref_edges, valid_src_edges, common_edges

        total_valid_ref_edges = len([n for n in ref_tree.traverse() if n.children and n.support > min_support_ref])
        result = {}
        if has_duplications:
            orig_target_size = len(source_tree)
            ntrees, ndups, sp_trees = source_tree.get_speciation_trees(
                autodetect_duplications=True, newick_only=True,
                target_attr=source_tree_attr, map_features=[source_tree_attr, "support"])

            if ntrees < max_treeko_splits_to_be_artifact:
                all_rf = []
                ref_found = []
                src_found = []
                tree_sizes = []
                all_max_rf = []
                common_names = 0

                for subtree_nw in sp_trees:

                    # if seedid and not use_collateral and (seedid not in subtree_nw):
                    #    continue
                    subtree = source_tree.__class__(subtree_nw, sp_naming_function=source_tree._speciesFunction)
                    if not subtree.children:
                        continue

                    # only necessary if rf function is going to filter by support
                    # value.  It slows downs the analysis, obviously, as it has to
                    # find the support for each node in the treeko tree from the
                    # original one.
                    if min_support_source > 0:
                        subtree_content = subtree.get_cached_content(store_attr='name')
                        for n in subtree.traverse():
                            if n.children:
                                n.support = source_tree.get_common_ancestor(subtree_content[n]).support

                    total_rf, max_rf, ncommon, valid_ref_edges, valid_src_edges, common_edges = _compare(subtree,
                                                                                                         ref_tree)

                    all_rf.append(total_rf)
                    all_max_rf.append(max_rf)
                    tree_sizes.append(ncommon)

                    if unrooted:
                        ref_found_in_src = len(common_edges) / float(len(valid_ref_edges)) if valid_ref_edges else None
                        src_found_in_ref = len(common_edges) / float(len(valid_src_edges)) if valid_src_edges else None
                    else:
                        # in rooted trees, we want to discount the root edge
                        # from the percentage of congruence. Otherwise we will never see a 0%
                        # congruence for totally different trees
                        ref_found_in_src = (len(common_edges) - 1) / float(len(valid_ref_edges) - 1) if len(
                            valid_ref_edges) > 1 else None
                        src_found_in_ref = (len(common_edges) - 1) / float(len(valid_src_edges) - 1) if len(
                            valid_src_edges) > 1 else None

                    if ref_found_in_src is not None:
                        ref_found.append(ref_found_in_src)

                    if src_found_in_ref is not None:
                        src_found.append(src_found_in_ref)

                if all_rf:
                    # Treeko speciation distance
                    alld = [_safe_div(all_rf[i], float(all_max_rf[i])) for i in range(len(all_rf))]
                    a = sum([alld[i] * tree_sizes[i] for i in range(len(all_rf))])
                    b = float(sum(tree_sizes))
                    treeko_d = a / b if a else 0.0
                    result["treeko_dist"] = treeko_d

                    result["rf"] = np.mean(all_rf)
                    result["max_rf"] = max(all_max_rf)
                    result["effective_tree_size"] = np.mean(tree_sizes)
                    result["norm_rf"] = np.mean(
                        [_safe_div(all_rf[i], float(all_max_rf[i])) for i in range(len(all_rf))])

                    result["ref_edges_in_source"] = np.mean(ref_found)
                    result["source_edges_in_ref"] = np.mean(src_found)

                    result["source_subtrees"] = len(all_rf)
                    result["common_edges"] = set()
                    result["source_edges"] = set()
                    result["ref_edges"] = set()
        else:
            total_rf, max_rf, ncommon, valid_ref_edges, valid_src_edges, common_edges = _compare(source_tree, ref_tree)

            result["rf"] = float(total_rf) if max_rf else "NA"
            result["max_rf"] = float(max_rf)
            if unrooted:
                result["ref_edges_in_source"] = len(common_edges) / float(
                    len(valid_ref_edges)) if valid_ref_edges else "NA"
                result["source_edges_in_ref"] = len(common_edges) / float(
                    len(valid_src_edges)) if valid_src_edges else "NA"
            else:
                # in rooted trees, we want to discount the root edge from the
                # percentage of congruence. Otherwise we will never see a 0%
                # congruence for totally different trees
                result["ref_edges_in_source"] = (len(common_edges) - 1) / float(len(valid_ref_edges) - 1) if len(
                    valid_ref_edges) > 1 else "NA"
                result["source_edges_in_ref"] = (len(common_edges) - 1) / float(len(valid_src_edges) - 1) if len(
                    valid_src_edges) > 1 else "NA"

            result["effective_tree_size"] = ncommon
            result["norm_rf"] = total_rf / float(max_rf) if max_rf else "NA"
            result["treeko_dist"] = "NA"
            result["source_subtrees"] = 1
            result["common_edges"] = common_edges
            result["source_edges"] = valid_src_edges
            result["ref_edges"] = valid_ref_edges
        return result

    def _diff(self, t2, output='topology', attr_t1='name', attr_t2='name', color=True):
        """
        .. versionadded:: 2.3

        Show or return the difference between two tree topologies.

        :param [raw|table|topology|diffs|diffs_tab] output: Output type

        """
        from ..tools import ete_diff
        difftable = ete_diff.treediff(self, t2, attr1=attr_t1, attr2=attr_t2)
        if output == "topology":
            ete_diff.show_difftable_topo(difftable, attr_t1, attr_t2, usecolor=color)
        elif output == "diffs":
            ete_diff.show_difftable(difftable)
        elif output == "diffs_tab":
            ete_diff.show_difftable_tab(difftable)
        elif output == 'table':
            rf, rf_max, _, _, _, _, _ = self.robinson_foulds(t2, attr_t1=attr_t1, attr_t2=attr_t2)[:2]
            ete_diff.show_difftable_summary(difftable, rf, rf_max)
        else:
            return difftable

    def iter_edges(self, cached_content=None):
        '''
        .. versionadded:: 2.3

        Iterate over the list of edges of a tree. Each edge is represented as a
        tuple of two elements, each containing the list of nodes separated by
        the edge.
        '''

        if not cached_content:
            cached_content = self.get_cached_content()
        all_leaves = cached_content[self]
        for n, side1 in cached_content.items():
            yield (side1, all_leaves - side1)

    def get_edges(self, cached_content=None):
        '''
        .. versionadded:: 2.3

        Returns the list of edges of a tree. Each edge is represented as a
        tuple of two elements, each containing the list of nodes separated by
        the edge.
        '''

        return [edge for edge in self.iter_edges(cached_content)]

    def standardize(self, delete_orphan=True, preserve_branch_length=True):
        """
        .. versionadded:: 2.3

        process current tree structure to produce a standardized topology: nodes
        with only one child are removed and multifurcations are automatically resolved.


        """
        self.resolve_polytomy()

        for n in self.get_descendants():
            if len(n.children) == 1:
                n.delete(prevent_nondicotomic=True,
                         preserve_branch_length=preserve_branch_length)

    def get_topology_id(self, attr="name"):
        '''
        .. versionadded:: 2.3

        Returns the unique ID representing the topology of the current tree. Two
        trees with the same topology will produce the same id. If trees are
        unrooted, make sure that the root node is not binary or use the
        tree.unroot() function before generating the topology id.

        This is useful to detect the number of unique topologies over a bunch of
        trees, without requiring full distance methods.

        The id is, by default, calculated based on the terminal node's names. Any
        other node attribute could be used instead.


        '''
        edge_keys = []
        for s1, s2 in self.get_edges():
            k1 = sorted([getattr(e, attr) for e in s1])
            k2 = sorted([getattr(e, attr) for e in s2])
            edge_keys.append(sorted([k1, k2]))
        return md5(str(sorted(edge_keys)).encode('utf-8')).hexdigest()

    # def get_partitions(self):
    #     """
    #     .. versionadded: 2.1

    #     It returns the set of all possible partitions under a
    #     node. Note that current implementation is quite inefficient
    #     when used in very large trees.

    #     t = Tree("((a, b), e);")
    #     partitions = t.get_partitions()

    #     # Will return:
    #     # a,b,e
    #     # a,e
    #     # b,e
    #     # a,b
    #     # e
    #     # b
    #     # a
    #     """
    #     all_leaves = frozenset(self.get_leaf_names())
    #     all_partitions = set([all_leaves])
    #     for n in self.iter_descendants():
    #         p1 = frozenset(n.get_leaf_names())
    #         p2 = frozenset(all_leaves - p1)
    #         all_partitions.add(p1)
    #         all_partitions.add(p2)
    #     return all_partitions

    def convert_to_ultrametric(self, tree_length=None, strategy='balanced'):
        """
        .. versionadded: 2.1

        Converts a tree into ultrametric topology (all leaves must have
        the same distance to root). Note that, for visual inspection
        of ultrametric trees, node.img_style["size"] should be set to
        0.
        """

        # Could something like this replace the old algorithm?
        # most_distant_leaf, tree_length = self.get_farthest_leaf()
        # for leaf in self:
        #    d = leaf.get_distance(self)
        #    leaf.dist += (tree_length - d)
        # return

        # pre-calculate how many splits remain under each node
        node2max_depth = {}
        for node in self.traverse("postorder"):
            if not node.is_leaf():
                max_depth = max([node2max_depth[c] for c in node.children]) + 1
                node2max_depth[node] = max_depth
            else:
                node2max_depth[node] = 1
        node2dist = {self: 0.0}
        if not tree_length:
            most_distant_leaf, tree_length = self.get_farthest_leaf()
        else:
            tree_length = float(tree_length)

        step = tree_length / node2max_depth[self]
        for node in self.iter_descendants("levelorder"):
            if strategy == "balanced":
                node.dist = (tree_length - node2dist[node.up]) / node2max_depth[node]
                node2dist[node] = node.dist + node2dist[node.up]
            elif strategy == "fixed":
                if not node.is_leaf():
                    node.dist = step
                else:
                    node.dist = tree_length - ((node2dist[node.up]) * step)
                node2dist[node] = node2dist[node.up] + 1
            node.dist = node.dist

    def check_monophyly(self, values, target_attr, ignore_missing=False,
                        unrooted=False):
        """
        .. versionadded: 2.2

        Returns True if a given target attribute is monophyletic under
        this node for the provided set of values.

        If not all values are represented in the current tree
        structure, a ValueError exception will be raised to warn that
        strict monophyly could never be reached (this behaviour can be
        avoided by enabling the `ignore_missing` flag.

        :param values: a set of values for which monophyly is
            expected.

        :param target_attr: node attribute being used to check
            monophyly (i.e. species for species trees, names for gene
            family trees, or any custom feature present in the tree).

        :param False ignore_missing: Avoid raising an Exception when
            missing attributes are found.


        .. versionchanged: 2.3

        :param False unrooted: If True, tree will be treated as unrooted, thus
          allowing to find monophyly even when current outgroup is splitting a
          monophyletic group.

        :returns: the following tuple
                  IsMonophyletic (boolean),
                  clade type ('monophyletic', 'paraphyletic' or 'polyphyletic'),
                  leaves breaking the monophyly (set)

        """

        if type(values) != set:
            values = set(values)

        # This is the only time I traverse the tree, then I use cached
        # leaf content
        n2leaves = self.get_cached_content()

        # Raise an error if requested attribute values are not even present
        if ignore_missing:
            found_values = set([getattr(n, target_attr) for n in n2leaves[self]])
            missing_values = values - found_values
            values = values & found_values

        # Locate leaves matching requested attribute values
        targets = set([leaf for leaf in n2leaves[self]
                       if getattr(leaf, target_attr) in values])
        if not ignore_missing:
            if values - set([getattr(leaf, target_attr) for leaf in targets]):
                raise ValueError(
                    'The monophyly of the provided values could never be reached, as not all of them exist in the tree.'
                    ' Please check your target attribute and values, or set the ignore_missing flag to True')

        if unrooted:
            smallest = None
            for side1, side2 in self.iter_edges(cached_content=n2leaves):
                if targets.issubset(side1) and (not smallest or len(side1) < len(smallest)):
                    smallest = side1
                elif targets.issubset(side2) and (not smallest or len(side2) < len(smallest)):
                    smallest = side2
                if smallest is not None and len(smallest) == len(targets):
                    break
            foreign_leaves = smallest - targets
        else:
            # Check monophyly with get_common_ancestor. Note that this
            # step does not require traversing the tree again because
            # targets are node instances instead of node names, and
            # get_common_ancestor function is smart enough to detect it
            # and avoid unnecessary traversing.
            common = self.get_common_ancestor(targets)
            observed = n2leaves[common]
            foreign_leaves = set([leaf for leaf in observed
                                  if getattr(leaf, target_attr) not in values])

        if not foreign_leaves:
            return True, "monophyletic", foreign_leaves
        else:
            # if the requested attribute is not monophyletic in this
            # node, let's differentiate between poly and paraphyly.
            poly_common = self.get_common_ancestor(foreign_leaves)
            # if the common ancestor of all foreign leaves is self
            # contained, we have a paraphyly. Otherwise, polyphyly.
            polyphyletic = [leaf for leaf in poly_common if
                            getattr(leaf, target_attr) in values]
            if polyphyletic:
                return False, "polyphyletic", foreign_leaves
            else:
                return False, "paraphyletic", foreign_leaves

    def get_monophyletic(self, values, target_attr):
        """
        .. versionadded:: 2.2

        Returns a list of nodes matching the provided monophyly
        criteria. For a node to be considered a match, all
        `target_attr` values within and node, and exclusively them,
        should be grouped.

        :param values: a set of values for which monophyly is
            expected.

        :param target_attr: node attribute being used to check
            monophyly (i.e. species for species trees, names for gene
            family trees).

        """

        if type(values) != set:
            values = set(values)

        n2values = self.get_cached_content(store_attr=target_attr)

        is_monophyletic = lambda node: n2values[node] == values
        for match in self.iter_leaves(is_leaf_fn=is_monophyletic):
            if is_monophyletic(match):
                yield match

    def expand_polytomies(self, map_attr="name", polytomy_size_limit=5,
                          skip_large_polytomies=False):
        '''
        .. versionadded:: 2.3

        Given a tree with one or more polytomies, this functions returns the
        list of all trees (in newick format) resulting from the combination of
        all possible solutions of the multifurcated nodes.

        .. warning:

           Please note that the number of of possible binary trees grows
           exponentially with the number and size of polytomies. Using this
           function with large multifurcations is not feasible:

           polytomy size: 3 number of binary trees: 3
           polytomy size: 4 number of binary trees: 15
           polytomy size: 5 number of binary trees: 105
           polytomy size: 6 number of binary trees: 945
           polytomy size: 7 number of binary trees: 10395
           polytomy size: 8 number of binary trees: 135135
           polytomy size: 9 number of binary trees: 2027025

        http://ajmonline.org/2010/darwin.php
        '''

        class TipTuple(tuple):
            pass

        def add_leaf(tree, label):
            yield (label, tree)
            if not isinstance(tree, TipTuple) and isinstance(tree, tuple):
                for left in add_leaf(tree[0], label):
                    yield (left, tree[1])
                for right in add_leaf(tree[1], label):
                    yield (tree[0], right)

        def enum_unordered(labels):
            if len(labels) == 1:
                yield labels[0]
            else:
                for tree in enum_unordered(labels[1:]):
                    for new_tree in add_leaf(tree, labels[0]):
                        yield new_tree

        n2subtrees = {}
        for n in self.traverse("postorder"):
            if n.is_leaf():
                subtrees = [getattr(n, map_attr)]
            else:
                subtrees = []
                if len(n.children) > polytomy_size_limit:
                    if skip_large_polytomies:
                        for childtrees in itertools.product(*[n2subtrees[ch] for ch in n.children]):
                            subtrees.append(TipTuple(childtrees))
                    else:
                        raise TreeError("Found polytomy larger than current limit: %s" % n)
                else:
                    for childtrees in itertools.product(*[n2subtrees[ch] for ch in n.children]):
                        subtrees.extend([TipTuple(subtree) for subtree in enum_unordered(childtrees)])

            n2subtrees[n] = subtrees
        return ["%s;" % str(nw) for nw in n2subtrees[self]]  # tuples are in newick format ^_^

    def resolve_polytomy(self, default_dist=0.0, default_support=0.0,
                         recursive=True):
        """
        .. versionadded: 2.2

        Resolve all polytomies under current node by creating an
        arbitrary dicotomic structure among the affected nodes. This
        function randomly modifies current tree topology and should
        only be used for compatibility reasons (i.e. programs
        rejecting multifurcated node in the newick representation).

        :param 0.0 default_dist: artificial branch distance of new
            nodes.

        :param 0.0 default_support: artificial branch support of new
            nodes.

        :param True recursive: Resolve any polytomy under this
             node. When False, only current node will be checked and fixed.
        """

        def _resolve(node):
            if len(node.children) > 2:
                children = list(node.children)
                node.children = []
                next_node = root = node
                for i in range(len(children) - 2):
                    next_node = next_node.add_child()
                    next_node.dist = default_dist
                    next_node.support = default_support

                next_node = root
                for ch in children:
                    next_node.add_child(ch)
                    if ch != children[-2]:
                        next_node = next_node.children[0]

        target = [self]
        if recursive:
            target.extend([n for n in self.get_descendants()])
        for n in target:
            _resolve(n)

    def add_face(self, face, column, position="branch-right"):
        """
        .. versionadded: 2.1

        Add a fixed face to the node.  This type of faces will be
        always attached to nodes, independently of the layout
        function.

        :argument face: a Face or inherited instance
        :argument column: An integer number starting from 0
        :argument "branch-right" position: Posible values are:
          "branch-right", "branch-top", "branch-bottom", "float",
          "aligned"
        """

        if not hasattr(self, "_faces"):
            self._faces = _FaceAreas()

        if position not in FACE_POSITIONS:
            raise ValueError("face position not in %s" % FACE_POSITIONS)

        if isinstance(face, Face):
            getattr(self._faces, position).add_face(face, column=column)
        else:
            raise ValueError("not a Face instance")


    @staticmethod
    def from_parent_child_table(parent_child_table):
        """Converts a parent-child table into an ETE Tree instance.

        :argument parent_child_table: a list of tuples containing parent-child
           relationships. For example: [("A", "B", 0.1), ("A", "C", 0.2), ("C",
           "D", 1), ("C", "E", 1.5)]. Where each tuple represents: [parent, child,
           child-parent-dist]

        :returns: A new Tree instance

        :example:

        >>> tree = Tree.from_parent_child_table([("A", "B", 0.1), ("A", "C", 0.2), ("C", "D", 1), ("C", "E", 1.5)])
        >>> print tree

        """

        def get_node(nodename, dist=None):
            if nodename not in nodes_by_name:
                nodes_by_name[nodename] = Tree(name=nodename, dist=dist)
            node = nodes_by_name[nodename]
            if dist is not None:
                node.dist = dist
            node.name = nodename
            return nodes_by_name[nodename]

        nodes_by_name = {}
        for columns in parent_child_table:
            if len(columns) == 3:
                parent_name, child_name, distance = columns
                dist = float(distance)
            else:
                parent_name, child_name = columns
                dist = None
            parent = get_node(parent_name)
            parent.add_child(get_node(child_name, dist=dist))

        root = parent.get_tree_root()
        return root

    @staticmethod
    def from_skbio(skbio_tree, map_attributes=None):
        """Converts a scikit-bio TreeNode object into ETE Tree object.

        :argument skbio_tree: a scikit bio TreeNode instance

        :argument None map_attributes: A list of attribute nanes in the
           scikit-bio tree that should be mapped into the ETE tree
           instance. (name, id and branch length are always mapped)

        :returns: A new Tree instance

        :example:

        >>> tree = Tree.from_skibio(skbioTree, map_attributes=["value"])

        """
        from skbio import TreeNode as skbioTreeNode

        def get_ete_node(skbio_node):
            ete_node = all_nodes.get(skbio_node, Tree())
            if skbio_node.length is not None:
                ete_node.dist = float(skbio_node.length)
            ete_node.name = skbio_node.name
            ete_node.add_features(id=skbio_node.id)
            if map_attributes:
                for a in map_attributes:
                    ete_node.add_feature(a, getattr(skbio_node, a, None))
            return ete_node

        all_nodes = {}
        if isinstance(skbio_tree, skbioTreeNode):
            for node in skbio_tree.preorder(include_self=True):
                all_nodes[node] = get_ete_node(node)
                ete_node = all_nodes[node]
                for ch in node.children:
                    ete_ch = get_ete_node(ch)
                    ete_node.add_child(ete_ch)
                    all_nodes[ch] = ete_ch
            return ete_ch.get_tree_root()

    def phonehome(self):
        from .. import _ph
        _ph.call()


def _translate_nodes(root, *nodes):
    name2node = dict([[n, None] for n in nodes if type(n) is str])
    for n in root.traverse():
        if n.name in name2node:
            if name2node[n.name] is not None:
                raise TreeError("Ambiguous node name: " + str(n.name))
            else:
                name2node[n.name] = n

    if None in list(name2node.values()):
        notfound = [key for key, value in name2node.items() if value is None]
        raise ValueError("Node names not found: " + str(notfound))

    valid_nodes = []
    for n in nodes:
        if type(n) is not str:
            if type(n) is not root.__class__:
                raise TreeError("Invalid target node: " + str(n))
            else:
                valid_nodes.append(n)

    valid_nodes.extend(list(name2node.values()))
    if len(valid_nodes) == 1:
        return valid_nodes[0]
    else:
        return valid_nodes


class HsTree(TreeNode):
#    USE_NODE_DICT = True

    def __init__(tree, *args, **kwa):
        if args:
            newick = args[0]
        else:
            try:
                newick = kwa['newick']
            except KeyError:
                newick = None

        super(HsTree, tree).__init__(*args, **kwa)
        try:
            tree.mass_migrations = kwa['mass_migrations']
        except KeyError:
            tree.mass_migrations = []
            if newick is not None:
                if os.path.exists(newick+'.mass_migrations.json'):
                    mass_migration_ls = json.load(open(newick+'.mass_migrations.json'))
                    for ms in mass_migration_ls:
                        tree.mass_migrations.append(MassMigration.from_dict(ms, tree))



#        if HsTree.USE_NODE_DICT:
 #           tree._create_node_dict()

#    def _create_node_dict(tree):
#        max_dist = max([tree.get_distance(l) for l in tree.get_leaves()])
#        tree.node_dict = {}
#        for node in tree.traverse():
#            node.time  = max_dist - tree.get_distance(node)
#            node.name = node.get_name()
#            tree.node_dict.update({node.name: node})

#    def get_node(tree, node_name):
#        return tree.node_dict[node_name]
    def get_time(node):
        rt = node.get_tree_root()
        max_dist = rt.get_farthest_leaf()[1]
        return max_dist - rt.get_distance(node)

    def get_name(tree):
        node1 = copy.deepcopy(tree)
        node1.sort_descendants()
        s = super(HsTree, node1).write(format=9)[:-1]
        #s = node1.write(format=9)[:-1]
        return s



    def write(tree, **kwa):
        """
        ete3 write method with the addition
        of writing mass migrations into
        separate json.

        :return: (ete3 write output, mass_migrations)

        """

        output = super().write(**kwa)
        mass_migration_ls = []
        for mm in tree.mass_migrations:
            md = mm.to_dict()
            mass_migration_ls.append(md)

        try:
            mm_filename = kwa['outfile'] + '.mass_migrations.json'
            with open(mm_filename, 'w') as f:
                json.dump(mass_migration_ls, f)
        except KeyError:
            pass

        return output, mass_migration_ls


    def add_mass_migration(tree, source, destination, fraction, time):
#        if HsTree.USE_NODE_DICT:
#            source = tree.node_dict[source_name]
#            destination = tree.node_dict[destination_name]
#        else:
#            for node in tree.traverse():
#                if node.get_name() == source_name:
#                    source = node
#                if node.get_name() == destination_name:
#                    destintion = node
        mm = MassMigration(source, destination, fraction, time)
        if mm not in tree.mass_migrations:
            tree.mass_migrations.append(mm)

    def add_property_to_nodes(tree, property_name, property_node_dict):
        """
        Adds the attribute property_name to nodes with newick
        given as property_node_dict keys and property values as
        dictionary values.

        TODO: Update for case if node_dict is available.

        Example:
        print tree

                 /-A1
              /-|
           /-|   \-A2
          |  |
        --|   \-B
          |
           \-C

        property='ne'
        property_node_dict={'(A2,A1);':4, 'C;':1}

        This adds the property node.ne to the nodes:
           /-A1
        --|
           \-A2
        --C
        """
        dic = {}
        for k,v in property_node_dict.iteritems():
            dic[newick_to_node_name(k)] = v

        for node in tree.traverse():
            node_name = node.get_name()
            try:
                setattr(node, property_name, dic[node_name])
            except KeyError:
                pass


    def add_properties_to_nodes(tree, properties, properties_node_dict):
        """
        Adds the attributes in the list properties to nodes with newick
        given as property_node_dict keys and a dictionary of
        {property:value} as dictionary values.

        Example:
        print tree

                 /-A1
              /-|
           /-|   \-A2
          |  |
        --|   \-B
          |
           \-C

        properties=['ne', 'color']
        property_node_dict={'(A2,A1);': {'ne':4, 'color': 'black'},
                             'C;': {'color': 'green'}}

        """
        dic = {}
        for k, v in properties_node_dict.iteritems():
            dic[newick_to_node_name(k)] = v


        for node in tree.traverse():
            node_name = node.get_name()
            for prop in properties:
                try:
                    setattr(node, prop, dic[node_name][prop])
                except KeyError:
                    pass

    def get_nodes_at_time(tree, time):
        nodes = []
        for n in tree.traverse():
            ancestors = n.get_ancestors()
            if ancestors:
                if ancestors[0].get_time() > time > n.get_time():
                    nodes.append(n)
        return nodes




    def plot(tree, ax=None, style='orthogonal',
                  node_name_fun=None,
                  node_name_format_fun=None,
                  leaf_name_fun=None,
                  leaf_name_format_fun=None,
                  line_format_fun=None,
                  migration_arrow_format_fun=None,
                  xtick_label_format_fun=None):
        """
        Plot ete tree.
        """
        default_node_format_args = dict(xycoords='data', ha='center',
                                    xytext=(0,1),
                                    textcoords='offset points',
                                    va='bottom',
                                    bbox=dict(boxstyle="round,pad=0.05", fc="w", alpha=0.5, lw=0),
                                        size=11)
        default_leaf_format_args = {'xytext':(5,0),
                                    'textcoords':'offset points',
                                    'va':'center'}
        default_line_format_args = {'color':'k'}
        default_migration_arrow_format_args = dict(arrowstyle="->, head_length = 0.5, head_width = .5",
                                                  color='r', linestyle='solid',linewidth=2,
                                                    zorder=-1)


        leaf_order = tree.get_leaf_names()

        if ax is None:
            fig = plt.figure(figsize=(12,len(leaf_order)*0.3))
            ax = plt.gca()

        assert style in ['orthogonal', 'diagonal']
        # don't plot node names if no function given
        if node_name_fun is None:
            node_name_fun = lambda node: False
        if node_name_format_fun is None:
            node_name_format_fun = lambda node: {}
        # plot leaf.name as leaf name by default
        if leaf_name_fun is None:
            leaf_name_fun = lambda node: node.name
        if leaf_name_format_fun is None:
            leaf_name_format_fun = lambda node: {}
        if line_format_fun is None:
            line_format_fun = lambda node: {}
        if migration_arrow_format_fun is None:
            migration_arrow_format_fun = lambda node: {}
        if xtick_label_format_fun is None:
            xtick_label_format_fun = lambda x, p: format(-int(x), ',')


        max_dist = max([tree.get_distance(l) for l in tree.get_leaves()])


        for i, node in enumerate(tree.traverse('postorder')):
            time =  node.get_time()
            if node.is_leaf():
                node.y = -leaf_order.index(node.name)
                leaf_name = leaf_name_fun(node)
                if leaf_name:
                    leaf_format_args = copy.deepcopy(default_leaf_format_args)
                    leaf_format_args.update(leaf_name_format_fun(node))
                    x = ax.annotate(leaf_name, xy=(-time, node.y),
                                xycoords='data', **leaf_format_args)

            else:
                l = node.children[0]
                r = node.children[1]
                node.y = (l.y+r.y)/2.

                for c in (l,r):

                    line_format_args = copy.deepcopy(default_line_format_args)
                    line_format_args.update(line_format_fun(c))
                    if style == 'orthogonal':
                        ax.hlines(c.y, -time, -c.get_time(), **line_format_args)
                        ax.vlines(-time,*sorted([c.y,node.y]), **line_format_args)

                    elif style == 'diagonal':
                        ax.plot([-time,-c.get_time()],[node.y, c.y])

                    if not c.is_leaf():
                        node_name = node_name_fun(c)
                        if node_name:
                            node_format_args = copy.deepcopy(default_node_format_args)
                            node_format_args.update(node_name_format_fun(c))
                            ax.annotate(node_name, xy=((-time-c.get_time())/2., c.y),
                                         **node_format_args)
                            #if "Haplochromis" in node_name:
                            #    print 'wtf'
                            #    return node

        for mm in tree.mass_migrations:
            #print "plotting migration one", mm.time, mm.source.get_name(), mm.destination.get_name()
            #ax.plot([-mm.time, -mm.time],sorted([mm.source.y, mm.destination.y]), color='r')
            #ax.arrow(-mm.time, mm.destination.y, 0 , mm.source.y - mm.destination.y,
            #                     length_includes_head=True, color='r', linestyle='dashed')
            migration_arrow_format_args = copy.deepcopy(default_migration_arrow_format_args)
            migration_arrow_format_args.update(migration_arrow_format_fun(c))
            ax.annotate("",xytext=(-mm.time, mm.destination.y), xy=(-mm.time,mm.source.y),
                         arrowprops=migration_arrow_format_args)
            ax.annotate("{}%".format(int(round(mm.fraction*100))), xy=(-mm.time, (mm.destination.y + mm.source.y)/2.),
                        ha='right',va='center', xytext=(-3,0),bbox=dict(boxstyle="round,pad=0.1", fc="w", alpha=0.5, lw=0),
                        textcoords='offset points', color='r')

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_yticks([])
        ax.set_ylabel('')
        ymin, ymax = ax.get_ylim()
        ax.set_ylim([ymin-(ymax-ymin)*0.05,ymax+(ymax-ymin)*0.01])
        ax.xaxis.tick_bottom()
        ax.get_xaxis().set_major_formatter(
                mpl.ticker.FuncFormatter(xtick_label_format_fun))
        return ax


    def plot_vertical(tree, ax=None, style='orthogonal',
                      orientation='left_to_right',
                      ax_pos=None,
                      origin=0,
                      use_distances=True,
                  node_name_fun=None,
                  node_name_format_fun=None,
                  leaf_name_fun=None,
                  leaf_name_format_fun=None,
                  line_format_fun=None,
                  migration_arrow_format_fun=None,
                  tick_label_format_fun=None):
        """
        Plot ete tree.
        """




        default_node_format_args = dict(xycoords='data', ha='center',
                                    xytext=(0,1),
                                    textcoords='offset points',
                                    va='bottom',
                                    bbox=dict(boxstyle="round,pad=0.05", fc="w", alpha=0.5, lw=0),
                                        size=11)
        default_leaf_format_args = {'textcoords':'offset points'}
                                    #, 'fontname':'monospace'
        default_line_format_args = {'color':'k'}
        default_migration_arrow_format_args = dict(arrowstyle="->, head_length = 0.5, head_width = .5",
                                                  color='r', linestyle='solid',linewidth=2,
                                                    zorder=-1)




        leaf_order = tree.get_leaf_names()

        if ax is None:
            fig = plt.figure(figsize=(12,len(leaf_order)*0.3))
            ax = plt.gca()





        assert style in ['orthogonal', 'diagonal']
        assert orientation in ['left_to_right', 'bottom_to_top',
                               'top_to_bottom', 'right_to_left']


        if orientation == 'left_to_right':
            specific_default_leaf_format = {'va':'center',
                                            'ha':'left'}
            edgeline = ax.hlines
            orthogonalline = ax.vlines
            sign = 1
            order = 1
            ax_pos1 = 'bottom'

        elif orientation == 'bottom_to_top':
            specific_default_leaf_format = {'va': 'bottom',
                                            'ha': 'center',
                                            'rotation': 90}
            edgeline = ax.vlines
            orthogonalline = ax.hlines
            sign = 1
            order = -1
            ax_pos1 = 'right'

        elif orientation == 'top_to_bottom':
            specific_default_leaf_format = {'va':'top',
                                        'ha': 'center',
                                            'rotation':90}
            edgeline = ax.vlines
            orthogonalline = ax.hlines
            sign = -1
            order = -1
            ax_pos1 = 'right'

        elif orientation == 'right_to_left':
            specific_default_leaf_format = {'va': 'center',
                                            'ha': 'right'}
            edgeline = ax.hlines
            orthogonalline = ax.vlines
            sign = -1
            order = 1
            ax_pos1 = 'bottom'


        if ax_pos is None:
            ax_pos = ax_pos1

        default_leaf_format_args.update(specific_default_leaf_format)

        default_leaf_format_args.update({'xytext': (sign * 5, 0)[::order]})


        # don't plot node names if no function given
        if node_name_fun is None:
            node_name_fun = lambda node: False
        if node_name_format_fun is None:
            node_name_format_fun = lambda node: {}
        # plot leaf.name as leaf name by default
        if leaf_name_fun is None:
            leaf_name_fun = lambda node: node.name
        if leaf_name_format_fun is None:
            leaf_name_format_fun = lambda node: {}
        if line_format_fun is None:
            line_format_fun = lambda node: {}
        if migration_arrow_format_fun is None:
            migration_arrow_format_fun = lambda node: {}
        if tick_label_format_fun is None:
            tick_label_format_fun = lambda x, p: format(-sign * int(x), ',')


        #max_dist = max([tree.get_distance(l) for l in tree.get_leaves()])


        for i, node in enumerate(tree.traverse('postorder')):
            time =  node.get_time()
            if node.is_leaf():
                node.y = origin -leaf_order.index(node.name)

                leaf_name = leaf_name_fun(node)
                if leaf_name:
                    leaf_format_args = copy.deepcopy(default_leaf_format_args)
                    leaf_format_args.update(leaf_name_format_fun(node))
                    x = ax.annotate(leaf_name, xy=(-sign * time , sign* node.y)[::order],
                                xycoords='data', **leaf_format_args)


            else:
                l = node.children[0]
                r = node.children[1]
                node.y =  (l.y+r.y)/2.

                for c in (l,r):

                    line_format_args = copy.deepcopy(default_line_format_args)
                    line_format_args.update(line_format_fun(c))
                    if style == 'orthogonal':

                        edgeline(sign * c.y, -sign * time, -sign * c.get_time(), **line_format_args)
                        orthogonalline(-sign * time, *sorted([sign * c.y, sign * node.y]), **line_format_args)

                    elif style == 'diagonal':
                        ax.plot([-time,-c.get_time()],[node.y, c.y])

                    if not c.is_leaf():
                        node_name = node_name_fun(c)
                        if node_name:
                            node_format_args = copy.deepcopy(default_node_format_args)
                            node_format_args.update(node_name_format_fun(c))
                            ax.annotate(node_name, xy=((-time-c.get_time())/2., c.y),
                                         **node_format_args)


        for mm in tree.mass_migrations:
            #print "plotting migration one", mm.time, mm.source.get_name(), mm.destination.get_name()
            #ax.plot([-mm.time, -mm.time],sorted([mm.source.y, mm.destination.y]), color='r')
            #ax.arrow(-mm.time, mm.destination.y, 0 , mm.source.y - mm.destination.y,
            #                     length_includes_head=True, color='r', linestyle='dashed')
            migration_arrow_format_args = copy.deepcopy(default_migration_arrow_format_args)
            migration_arrow_format_args.update(migration_arrow_format_fun(c))
            ax.annotate("",xytext=(-sign * mm.time, sign * mm.destination.y)[::order],
                        xy=(-sign * mm.time, sign * mm.source.y)[::order],
                         arrowprops=migration_arrow_format_args)

            ax.annotate("{}%".format(int(round(mm.fraction*100))),
                        xy=(-sign * mm.time, sign * (mm.destination.y + mm.source.y)/2.)[::order],
                         xytext=(-sign * 5,0)[::order],#ha='right',va='center',
                        bbox=dict(boxstyle="round,pad=0.1", fc="w", alpha=0.5, lw=0),
                        textcoords='offset points', color='r')



        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)




        if ax_pos:
            ax.spines[ax_pos].set_visible(True)
        else:
            ax.set_yticks([])
            ax.set_xticks([])


        if ax_pos in ['top','bottom']:
            ymin, ymax = ax.get_ylim()
            ax.set_ylim([ymin - (ymax - ymin) * 0.05, ymax + (ymax - ymin) * 0.01])
            #ax.xaxis.tick_bottom()
            ax.set_yticks([])
            ax.get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(tick_label_format_fun))
        elif ax_pos in ['left', 'right']:
            ymin, ymax = ax.get_xlim()
            ax.set_xlim([ymin - (ymax - ymin) * 0.05, ymax + (ymax - ymin) * 0.01])
            ax.yaxis.tick_right()
            ax.set_xticks([])
            ax.get_yaxis().set_major_formatter(
             mpl.ticker.FuncFormatter(tick_label_format_fun))



        #ymin, ymax = ax.get_ylim()
        #ax.set_ylim([ymin-(ymax-ymin)*0.05,ymax+(ymax-ymin)*0.01])
        #ax.xaxis.tick_bottom()
        #ax.get_xaxis().set_major_formatter(
        #        mpl.ticker.FuncFormatter(xtick_label_format_fun))
        #ax.set_xlim([1,10])


        return ax

    def plot_node_tree(tree, ax=None, style='orthogonal',
                      orientation='left_to_right',
                      ax_pos=None,
                      origin=0,
                      use_distances=True,
                      internal_node_to_present=False,
                  node_name_fun=None,
                  node_name_format_fun=None,
                  leaf_name_fun=None,
                  leaf_name_format_fun=None,
                  line_format_fun=None,
                  migration_arrow_format_fun=None,
                  tick_label_format_fun=None):
        """
        Plot ete tree.
        """




        default_node_format_args = dict(xycoords='data', ha='center',
                                    xytext=(0,1),
                                    textcoords='offset points',
                                    va='bottom',
                                    bbox=dict(boxstyle="round,pad=0.05", fc="w", alpha=0.5, lw=0),
                                        size=11)
        default_leaf_format_args = {'textcoords':'offset points'}
                                    #, 'fontname':'monospace'
        default_line_format_args = {'color':'k'}
        default_migration_arrow_format_args = dict(arrowstyle="->, head_length = 0.5, head_width = .5",
                                                  color='r', linestyle='solid',linewidth=2,
                                                    zorder=-1)



        if internal_node_to_present:
            leaf_order = [n.get_name() for n in tree.traverse('preorder')]
        else:
            leaf_order = tree.get_leaf_names()

        if ax is None:
            fig = plt.figure(figsize=(12, len(leaf_order)*0.3))
            ax = plt.gca()





        assert style in ['orthogonal', 'diagonal']
        assert orientation in ['left_to_right', 'bottom_to_top',
                               'top_to_bottom', 'right_to_left']


        if orientation == 'left_to_right':
            specific_default_leaf_format = {'va':'center',
                                            'ha':'left'}
            edgeline = ax.hlines
            orthogonalline = ax.vlines
            sign =  1
            order = 1
            ax_pos1 = 'bottom'


        elif orientation == 'bottom_to_top':
            specific_default_leaf_format = {'va': 'bottom',
                                            'ha': 'center',
                                            'rotation': 90}
            edgeline = ax.vlines
            orthogonalline = ax.hlines
            sign = 1
            order = -1
            ax_pos1 = 'right'


        elif orientation == 'top_to_bottom':
            specific_default_leaf_format = {'va':'top',
                                        'ha': 'center',
                                            'rotation': 90}
            edgeline = ax.vlines
            orthogonalline = ax.hlines
            sign = -1
            order = -1
            ax_pos1 = 'right'


        elif orientation == 'right_to_left':
            specific_default_leaf_format = {'va': 'center',
                                            'ha': 'right'}
            edgeline = ax.hlines
            orthogonalline = ax.vlines
            sign = -1
            order = 1
            ax_pos1 = 'bottom'
        else:
            raise ValueError("Orientation must be left_to_right,"
                             "bottom_to_top, top_to_bottom or right_to_left.")


        if ax_pos is None:
            ax_pos = ax_pos1

        default_leaf_format_args.update(specific_default_leaf_format)

        default_leaf_format_args.update({'xytext': (sign * 5, 0)[::order]})




        # don't plot node names if no function given
        if node_name_fun is None:
            node_name_fun = lambda node: False
        if node_name_format_fun is None:
            node_name_format_fun = lambda node: {}
        # plot leaf.name as leaf name by default
        if leaf_name_fun is None:
            leaf_name_fun = lambda node: node.name
        if leaf_name_format_fun is None:
            leaf_name_format_fun = lambda node: {}
        if line_format_fun is None:
            line_format_fun = lambda node: {}
        if migration_arrow_format_fun is None:
            migration_arrow_format_fun = lambda node: {}
        if tick_label_format_fun is None:
            tick_label_format_fun = lambda x, p: format(-sign * int(x), ',')


        #max_dist = max([tree.get_distance(l) for l in tree.get_leaves()])


        max_label_width = 0
        leaf_annots = []

        for i, node in enumerate(tree.traverse('postorder')):
            time =  node.get_time()
            if node.is_leaf():
                node.y = origin - leaf_order.index(node.name)

                leaf_name = leaf_name_fun(node)
                if leaf_name:
                    leaf_format_args = copy.deepcopy(default_leaf_format_args)
                    leaf_format_args.update(leaf_name_format_fun(node))
                    x = ax.annotate(leaf_name, xy=(-sign * time , sign* node.y)[::order],
                                xycoords='data', **leaf_format_args)
                    leaf_annots.append(x)

                    #print('y', bb.y0, bb.y1)
                    #print(inv.transform((bb.x1, bb.y1)))
                    #print(ax.transData.transform((bb.x1, bb.y1)))
                    #break

                    #print((bb.width,bb.height))
                    #print(ax.transData.transform((bb.width,bb.height)))




            else:
                l = node.children[0]
                r = node.children[1]

                if not internal_node_to_present:
                    node.y = (l.y+r.y)/2.
                else:
                    node.y = origin - leaf_order.index(node.get_name())


                #print(node.y)

                for c in (l,r):
                    line_format_args = copy.deepcopy(default_line_format_args)
                    line_format_args.update(line_format_fun(c))
                    if style == 'orthogonal':

                        edgeline(sign * (c.y), -sign * time, -sign * c.get_time(), **line_format_args)
                        orthogonalline(-sign * time, *sorted([sign * c.y , sign * node.y]), **line_format_args)

                    elif style == 'diagonal':
                        ax.plot([-time,-c.get_time()],[node.y, c.y])


                    if not c.is_leaf():
                        node_name = node_name_fun(c)
                        if node_name:
                            node_format_args = copy.deepcopy(default_node_format_args)
                            node_format_args.update(node_name_format_fun(c))
                            ax.annotate(node_name, xy=((-time-c.get_time())/2., c.y),
                                         **node_format_args)

        fig = plt.gcf()
        renderer = fig.canvas.get_renderer()
        inv = ax.transData.inverted()

        labeledges = []

        for x in leaf_annots:
            bb = x.get_window_extent(renderer=renderer)
            ((x0, y0), (x1, y1)) = inv.transform(bb)
            labeledges.append(((x0, y0), (x1, y1)) )


        labeledges = np.array(labeledges)
        #print(np.sort(labeledges[:, int((sign + 1) / 2), int((order - 1) / -2)])[::-order])


        edgecoord = np.sort(labeledges[:, int((sign + 1) / 2), int((order - 1) / -2)])[::-sign][0]
        #print(edgecoord)

        if internal_node_to_present:
            for node in tree.traverse():
                if not node.is_leaf():
                    edgeline(sign * node.y, -sign * node.get_time(),
                         edgecoord,
                         linestyle='dotted')#dotted


        for mm in tree.mass_migrations:
            #print "plotting migration one", mm.time, mm.source.get_name(), mm.destination.get_name()
            #ax.plot([-mm.time, -mm.time],sorted([mm.source.y, mm.destination.y]), color='r')
            #ax.arrow(-mm.time, mm.destination.y, 0 , mm.source.y - mm.destination.y,
            #                     length_includes_head=True, color='r', linestyle='dashed')
            migration_arrow_format_args = copy.deepcopy(default_migration_arrow_format_args)
            migration_arrow_format_args.update(migration_arrow_format_fun(c))
            ax.annotate("",xytext=(-sign * mm.time, sign * mm.destination.y)[::order],
                        xy=(-sign * mm.time, sign * mm.source.y)[::order],
                         arrowprops=migration_arrow_format_args)

            ax.annotate("{}%".format(int(round(mm.fraction*100))),
                        xy=(-sign * mm.time, sign * (mm.destination.y + mm.source.y)/2.)[::order],
                         xytext=(-sign * 5,0)[::order],#ha='right',va='center',
                        bbox=dict(boxstyle="round,pad=0.1", fc="w", alpha=0.5, lw=0),
                        textcoords='offset points', color='r')



        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)




        if ax_pos:
            ax.spines[ax_pos].set_visible(True)
        else:
            ax.set_yticks([])
            ax.set_xticks([])
            plt.tick_params(
                #axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                left=False,  # ticks along the bottom edge are off
                right=False,  # ticks along the top edge are off
                bottom=False,
                top=False,
                labelbottom=False)  # labels along the bottom edge are off



        ymin, ymax = ax.get_ylim()
        xmin, xmax = ax.get_xlim()


        if ax_pos in ['top','bottom']:
            ax.set_ylim([ymin - (ymax - ymin) * 0.05, ymax + (ymax - ymin) * 0.01])
            #ax.xaxis.tick_bottom()
            ax.set_yticks([])
            ax.get_xaxis().set_major_formatter(
                        mpl.ticker.FuncFormatter(tick_label_format_fun))

        elif ax_pos in ['left', 'right']:

            ax.set_xlim([ymin - (ymax - ymin) * 0.05, ymax + (ymax - ymin) * 0.01])
            ax.yaxis.tick_right()
            ax.set_xticks([])
            ax.get_yaxis().set_major_formatter(
                        mpl.ticker.FuncFormatter(tick_label_format_fun))



        if orientation == 'left_to_right':
            ax.set_xlim(xmin, max(xmax, edgecoord ))
        elif orientation == 'right_to_left':
            ax.set_xlim(min(xmin, edgecoord),xmax)
        elif orientation == 'top_to_bottom':
            ax.set_ylim(min(ymin, edgecoord), ymax)
        elif orientation == 'bottom_to_top':
            ax.set_ylim(ymin, max(ymax, edgecoord ))


        #ymin, ymax = ax.get_ylim()
        #ax.set_ylim([ymin-(ymax-ymin)*0.05,ymax+(ymax-ymin)*0.01])
        #ax.xaxis.tick_bottom()
        #ax.get_xaxis().set_major_formatter(
        #        mpl.ticker.FuncFormatter(xtick_label_format_fun))
        #ax.set_xlim([1,10])


        return ax



    def search_node_by_newick(tree, newick):
        for n in tree.traverse():
            if n.get_name() == newick:
                return n
        raise Exception('Node not found.')

    def set_leaf_order(tree, order, check_consistent=True):
        """
        Changes the tree so that the leaves
        conform to order.
        The order must be consistent with
        the branching structure.

        Parameters:
        tree ... ete3 tree object
        order ... list of leaf names

        Returns:
        None (tree changed in place)
        """
        order = list(order)
        for i, node in enumerate(tree.traverse('postorder')):
            if not node.is_leaf():
                l = node.children[0]
                r = node.children[1]
                lnames = l.get_leaf_names()
                rnames = r.get_leaf_names()
                if order.index(lnames[0]) > order.index(rnames[0]):
                    node.swap_children()
        if check_consistent:
            if tree.get_leaf_names() != order:
                raise Exception("The provided order is not conistent with tree: \nResulting order: {}\nInput order:{}".format(
                    tree.get_leaf_names(), order))


    def set_outgroup(tree, outgroup, end_at_present=True):
        """
        Set root of tree.

        :param outgroup: name of the outgroup leaf
        :param end_at_present: root at a point of the branch so that outgroup ends at time 0.
        :return: None
        """
        super().set_outgroup(outgroup)
        if end_at_present:
            outgroup_node = tree.search_nodes(name=outgroup)[0]
            ingroup_root = [n for n in tree.get_children() if n is not outgroup_node][0]
            time = outgroup_node.get_time()
            outgroup_node.dist = outgroup_node.dist + time / 2
            ingroup_root.dist = ingroup_root.dist - time / 2
            assert ingroup_root.dist > 0, \
                    "Outgroup branch too short to lead to present. Use end_at_present=False."

    def reverse(tree):
        tree.set_leaf_order(tree, tree.get_leaf_names()[::-1])




def align_fbranch_with_tree(fbranch, tree, outgroup):
    tree_no = copy.deepcopy(tree)
    # tree_no.ladderize()
    # remove outgroup
    tree_no.prune([n for n in tree_no.get_leaf_names() if n != outgroup])
    #
    fb = fbranch.copy()
    fb.index = fb.index.droplevel(0)
    fb = fb.drop(outgroup, axis=1)
    fb.index = [tuple(sorted(i.split(','))) for i in fb.index]
    row_order = []
    col_order = tree_no.get_leaf_names()
    for n in tree_no.iter_descendants(strategy='preorder'):
        row_order.append(tuple(sorted(n.get_leaf_names())))

    assert fb.shape == (len(row_order), len(col_order)), \
        "{} != ({},{}))".format(fb.shape, len(row_order), len(col_order))

    assert set(fb.index.values) == set(row_order), \
        "Samples in input tree and fbranch matrix not consistent: {} != {}".format(fb.index.values, row_order)

    assert set(fb.columns.values) == set(col_order), \
        "Samples in input tree and fbranch matrix not consistent: {} != {}".format(fb.columns.values, col_order)

    # order fb in a way that is conistent with the tree
    fb = fb.loc[row_order, col_order].iloc[::-1]  # ,::-1

    return fb, tree_no


def plot_fbranch(fbranch, tree_no_outgroup):
    n_rows, n_cols = fbranch.shape

    plt.rcParams['font.size'] = 12

    # depth = tree_sd2_1_ete.get_farthest_leaf(topology_only=True)[1] +2

    # visited = []

    # fig = plt.figure(figsize=(18,20))

    fig = plt.figure(figsize=(n_cols, n_rows * 0.5))

    ax0 = plt.subplot2grid((5, 7), (0, 2), rowspan=1, colspan=4)

    tree_no_outgroup.plot_node_tree(origin=-0.5, ax=ax0, orientation='top_to_bottom')
    ax0.axis('off')

    ax = plt.subplot2grid((5, 7), (1, 0), rowspan=4, colspan=2)

    tree_no_outgroup.plot_node_tree(origin=fbranch.shape[0] + 0.5, ax=ax, orientation='left_to_right',
                                    internal_node_to_present=True)
    ax.axis('off')
    # ax = plot_node_tree(tree_no_outgroup, ax=ax, x0=0,y0=len(branch_mat)-0.5,em=0.5,fontsize=12)

    # ax.set_xlabel('excess allele sharing with')
    # ax0.xaxis.set_label_position('top')

    # ax.set_xlim([-1,10])

    ax2 = plt.subplot2grid((5, 7), (1, 2), rowspan=4, colspan=4, sharey=ax, sharex=ax0)

    branch_mat0 = fbranch.copy()
    branch_mat0_masked = np.ma.array(branch_mat0, mask=np.isnan(branch_mat0))

    fmax = branch_mat0.max().max()
    fmin = branch_mat0.min().min()
    colors = np.concatenate([[[1, 1, 1, 1]], plt.cm.Reds(np.linspace(0., 1 * fmax / 0.15, 256))])
    mymap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colors)

    plt.pcolormesh(branch_mat0_masked, cmap=mymap, rasterized=True)  # ,cmap=jet

    # plt.scatter(zs.index.droplevel(0).values+0.5,zs.index.droplevel(1)+0.43,marker='*',s=10)

    ax2.set_xticks(np.arange(0.5, fbranch.shape[1]))
    # ax2.set_xticklabels(fbranch.columns,rotation=90)
    ax2.set_xticklabels([])
    ax2.set_yticks(np.arange(0.5, fbranch.shape[0]))
    ax2.set_yticklabels([])
    ax2.set_facecolor((0.85, 0.85, 0.85))
    plt.tight_layout()
    ax2.set_ylim(0, len(fbranch))
    plt.subplots_adjust(wspace=-0.1)
    plt.subplots_adjust(hspace=0.05)
    ax2.xaxis.tick_top()
    ax2.set_xlim([0, fbranch.shape[1]])

    for b in range(fbranch.shape[0]):
        l = ax2.axhline(y=b + 0.04, xmin=0, xmax=1, linewidth=1, color='grey', alpha=0.5)

    for b in range(fbranch.shape[1]):
        l = ax2.axvline(x=b - 0.02, ymin=0, ymax=1, linewidth=1, color='grey', alpha=0.5)

    cbaxes = fig.add_axes([ax2.get_position().xmax + 0.02, ax2.get_position().ymin,
                           0.03, ax2.get_position().ymax - ax2.get_position().ymin])
    # plt.subplot2grid((5,7), (1,6), rowspan=4,colspan=1)
    mappable = mpl.cm.ScalarMappable(cmap=mymap)
    mappable.set_array([fmin, fmax])
    cbar = plt.colorbar(mappable, cax=cbaxes, label='$f_b$')
    cbar.set_label('$f_b$', size=24)
    cbar.solids.set_rasterized(True)

    plt.axis('tight')

    plt.subplots_adjust(wspace=0.05)

    plt.tight_layout()

    return fig.gca()

def main():
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        description="Plot f-branch statistic as produced by Dsuite. "
                                                    "Produces .png and .svg files.",
                                                                                 add_help=True)
    argparser.add_argument("fbranch", type=argparse.FileType('r'),
                            help="Path to file containing f-branch matrix as produced by Dsuite Fbranch.")
    argparser.add_argument("tree", type = argparse.FileType('r'),
                            help="Path to .newick tree file as given to Dsuite Fbranch.")

    argparser.add_argument("-n", "--run-name", type=str,
                            help="Base file name for output plots.",default="fbranch")
    argparser.add_argument( "--outgroup", type=str,
                        help="Outgroup name in newick file.",default="Outgroup")

    args = argparser.parse_args()
    fb = pd.read_csv(args.fbranch, sep='\t',
                  index_col=[0,1])
    tree = HsTree(args.tree.read())
    fb1, tree_no_outgroup = align_fbranch_with_tree(fb, tree,
                                                     outgroup=args.outgroup)
    ax = plot_fbranch(fb1, tree_no_outgroup)
    plt.savefig(args.run_name+'.svg', bbox_inches='tight')
    plt.savefig(args.run_name+'.png', bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    sys.exit(main())
