# -*- coding: utf-8 -*-
"""
    sphinx.domains.chm
    ~~~~~~~~~~~~~~~~~~~

    The CHARMM domain.

    :copyright: This is a work by the US Government, and is not subject to copyright.

    Defines directives for sphinx domain:
        keyword - kw
        command - cmd
        developer - dev
        spec - spec
"""

from docutils import nodes
from docutils.parsers.rst import directives

from sphinx import addnodes
from sphinx.roles import XRefRole
from sphinx.locale import l_, _
from sphinx.domains import Domain, ObjType, Index
from sphinx.directives import ObjectDescription
from sphinx.util.nodes import make_refnode
from sphinx.util.compat import Directive
from sphinx.util.docfields import Field, GroupedField, TypedField, DocFieldTransformer


def CHARmmify(string, len=4):
    string = string.strip()
    pre = string[:len].upper()
    post = string[len:].lower()
    return pre + post


class ChmObject(ObjectDescription):
    has_content =True
    required_arguments = 1
    optional_arguments = 10
    final_argument_whitespace = False

    def add_target_and_index(self, name, sig, signode):
        targetname = self.objtype + '-' + name
        if targetname not in self.state.document.ids:
            signode['names'].append(targetname)
            signode['ids'].append(targetname)
            signode['first'] = (not self.names)
            self.state.document.note_explicit_target(signode)

            objects = self.env.domaindata['chm']['objects']
            key = (self.objtype, name)
            if key in objects:
                self.state_machine.reporter.warning(
                    'duplicate description of %s %s, ' % (self.objtype, name) +
                    'other instance in ' + self.env.doc2path(objects[key]),
                    line=self.lineno)
            objects[key] = self.env.docname
        indextext = self.get_index_text(self.objtype, name)
        if indextext:
            self.indexnode['entries'].append(('single', indextext,
                                              targetname, ''))

    def get_index_text(self, objectname, name):
        if self.objtype == 'keyword':
            return _('%s (keyword)') % name
        elif self.objtype == 'command':
            return _('%s (command)') % name
        elif self.objtype == 'developer':
            return _('%s (developer)') % name
        elif self.objtype == 'spec':
            return _('%s (spec)') % name
        return ''


class KwDirective(ChmObject):
    option_spec = {
        'noindex': directives.flag,
        'default': directives.unchanged,
        'length': directives.positive_int,
    }

    def handle_signature(self, sig, signode):
        name = CHARmmify(sig, self.length) + ' '
        signode += addnodes.desc_name(name, name)
        for arg in self.arguments[1:]:
            arg = u' %s ' % arg
            signode += addnodes.desc_type(arg, arg)
        return sig.lower()[:self.length]

    def run(self):
        """
        Main directive entry function, called by docutils upon encountering the
        directive.

        This directive is meant to be quite easily subclassable, so it delegates
        to several additional methods.  What it does:

        * find out if called as a domain-specific directive, set self.domain
        * create a `desc` node to fit all description inside
        * parse standard options, currently `noindex`
        * create an index node if needed as self.indexnode
        * parse all given signatures (as returned by self.get_signatures())
          using self.handle_signature(), which should either return a name
          or raise ValueError
        * add index entries using self.add_target_and_index()
        * parse the content and handle doc fields in it
        """
        if ':' in self.name:
            self.domain, self.objtype = self.name.split(':', 1)
        else:
            self.domain, self.objtype = '', self.name
        self.env = self.state.document.settings.env
        self.indexnode = addnodes.index(entries=[])

        self.length = self.options.get('length', 4)
        node = addnodes.desc()
        node.document = self.state.document
        node['domain'] = self.domain
        # 'desctype' is a backwards compatible attribute
        node['objtype'] = node['desctype'] = self.objtype
        node['noindex'] = noindex = ('noindex' in self.options)

        self.names = []
        signatures = self.get_signatures()
        for i, sig in enumerate(signatures):
            # add a signature node for each signature in the current unit
            # and add a reference target for it
            signode = addnodes.desc_signature(sig, '')
            signode['first'] = False
            node.append(signode)
            try:
                # name can also be a tuple, e.g. (classname, objname);
                # this is strictly domain-specific (i.e. no assumptions may
                # be made in this base class)
                name = self.handle_signature(sig, signode)
                #print name
            except ValueError:
                # signature parsing failed
                signode.clear()
                signode += addnodes.desc_name(sig, sig)
                continue  # we don't want an index entry here
            if name not in self.names:
                self.names.append(name)
                if not noindex:
                    # only add target and index entry if this is the first
                    # description of the object with this name in this desc block
                    self.add_target_and_index(name, sig, signode)

        #####################
        # options processing
        #####################

        # default
        try:
            mydef = self.options['default']
            valnode = nodes.inline(mydef, mydef)
            boldnode = nodes.strong('Default: ', 'Default: ')
            #defnode = nodes.paragraph()
            defnode = addnodes.desc_content()
            defnode += boldnode
            defnode += valnode
            node += defnode
        except KeyError:
            pass

        #####################
        # content processing
        #####################

        contentnode = addnodes.desc_content()
        node.append(contentnode)
        # build node
        if self.names:
            # needed for association of version{added,changed} directives
            self.env.temp_data['object'] = self.names[0]
        self.before_content()
        self.state.nested_parse(self.content, self.content_offset, contentnode)
        DocFieldTransformer(self).transform_all(contentnode)
        self.env.temp_data['object'] = None
        self.after_content()
        return [self.indexnode, node]


class CmdDirective(KwDirective):
    option_spec = {
        'noindex': directives.flag,
        'default': directives.unchanged,
    }


class DevDirective(ChmObject):
    pass
    #final_argument_whitespace = False

    #option_spec = {
    #    'noindex': directives.flag,
    #    'default': directives.unchanged,
    #}

    #def handle_signature(self, sig, signode):
    #    ssig = sig.split()
    #    name = u'\uFEFF     %s-spec    \uFEFF' % ssig[0]
    #    signode += addnodes.desc_name(name, name)
    #    for arg in self.arguments[1:]:
    #        arg = ' ' + arg
    #        signode += addnodes.desc_type(arg, arg)
    #    return sig.lower()

class SpecDirective(ChmObject):

    def handle_signature(self, sig, signode):
        ssig = sig.split()
        name = u'%s-spec ' % ssig[0]
        signode += addnodes.desc_name(name, name)
        for arg in self.arguments[1:]:
            arg = u' %s ' % arg
            signode += addnodes.desc_type(arg, arg)
        return sig.lower()


class ChmXRefRole(XRefRole):
    def process_link(self, env, refnode, has_explicit_title, title, target):
        if not has_explicit_title:
            title = title.lstrip('.')   # only has a meaning for the target
            target = target.lstrip('~') # only has a meaning for the title
            # if the first character is a tilde, don't display the module/class
            # parts of the contents
            if title[0:1] == '~':
                title = title[1:]
                dot = title.rfind('.')
                if dot != -1:
                    title = title[dot+1:]
        # if the first character is a dot, search more specific namespaces first
        # else search builtins first
        if target[0:1] == '.':
            target = target[1:]
            refnode['refspecific'] = True
        #print title, type(title)
        #print target, type(target)
        ##print refnode.__dict__
        ##print env.domaindata['chm']['objects']
        #print env.domains['chm'].__dict__
        #print env.app.domains['chm'].__dict__
        #print env.todo_all_todos[0]
        return title, target


class DevXRefRole(ChmXRefRole):
    def process_link(self, env, refnode, has_explicit_title, title, target):
        refnode['chm:developer'] = env.temp_data.get('chm:developer')
        title, target = super(DevXRefRole, self).process_link(env, refnode, has_explicit_title, title, target)
        title = CHARmmify(title)
        return title, target


class KwXRefRole(ChmXRefRole):
    def process_link(self, env, refnode, has_explicit_title, title, target):
        refnode['chm:keyword'] = env.temp_data.get('chm:keyword')
        target = target.lower()
        title, target = super(KwXRefRole, self).process_link(env, refnode, has_explicit_title, title, target)
        title = CHARmmify(title)
        return title, target


class CmdXRefRole(ChmXRefRole):
    def process_link(self, env, refnode, has_explicit_title, title, target):
        refnode['chm:command'] = env.temp_data.get('chm:command')
        target = target.lower()
        title, target = super(CmdXRefRole, self).process_link(env, refnode, has_explicit_title, title, target)
        title = CHARmmify(title)
        return title, target


class SpecXRefRole(ChmXRefRole):
    def process_link(self, env, refnode, has_explicit_title, title, target):
        refnode['chm:spec'] = env.temp_data.get('chm:spec')
        target = target.lower()
        #title, target = super(SpecXRefRole, self).process_link(env, refnode, has_explicit_title, title, target)
        title, target = super(SpecXRefRole, self).process_link(env, refnode, True, title, target)
        title = '%s-spec' % title
        return title, target


class ChmKwIndex(Index):
    name = 'kwindex'
    localname = l_('CHARMM Keywords Index')
    shortname = l_('keywords')

    def generate(self, docnames=None):
        content = {}
        # list of prefixes to ignore
        ignores = self.domain.env.config['kwindex_common_prefix']
        ignores = sorted(ignores, key=len, reverse=True)
        # list of all keywordss, sorted by keyword name
        keywords = sorted(self.domain.data['keywords'].iteritems(),
                         key=lambda x: x[0].lower())
        # sort out collapsable keywords
        prev_kwname = ''
        num_toplevels = 0
        for kwname, (docname, synopsis, platforms, deprecated) in keywords:
            if docnames and docname not in docnames:
                continue

            for ignore in ignores:
                if kwname.startswith(ignore):
                    kwname = kwname[len(ignore):]
                    stripped = ignore
                    break
            else:
                stripped = ''

            # we stripped the whole keyword name?
            if not kwname:
                kwname, stripped = stripped, ''

            entries = content.setdefault(kwname[0].lower(), [])

            package = kwname.split('.')[0]
            if package != kwname:
                # it's a subkeyword
                if prev_kwname == package:
                    # first subkeyword - make parent a group head
                    if entries:
                        entries[-1][1] = 1
                elif not prev_kwname.startswith(package):
                    # subkeyword without parent in list, add dummy entry
                    entries.append([stripped + package, 1, '', '', '', '', ''])
                subtype = 2
            else:
                num_toplevels += 1
                subtype = 0

            qualifier = deprecated and _('Deprecated') or ''
            entries.append([stripped + kwname, subtype, docname,
                            'keyword-' + stripped + kwname, platforms,
                            qualifier, synopsis])
            prev_kwname = kwname

        # apply heuristics when to collapse modindex at page load:
        # only collapse if number of toplevel keywords is larger than
        # number of subkeywords
        collapse = len(keywords) - num_toplevels < num_toplevels

        # sort by first letter
        content = sorted(content.iteritems())

        return content, collapse


class ChmDomain(Domain):
    name = 'chm'
    label = 'CHARMM'
    object_types = {
        'keyword': ObjType(l_('keyword'), 'kw', 'obj'),
        'command': ObjType(l_('command'), 'cmd', 'obj'),
        'developer': ObjType(l_('developer'), 'dev', 'obj'),
        'spec': ObjType(l_('spec'), 'spec', 'obj'),
    }

    directives = {
        'keyword': KwDirective,
        'command': CmdDirective,
        'developer': DevDirective,
        'spec': SpecDirective,
    }

    roles = {
        'kw': KwXRefRole(),
        'cmd': CmdXRefRole(),
        'dev': ChmXRefRole(),
        'spec': SpecXRefRole(),
    }

    initial_data = {
        'objects': {},  # fullname -> docname, objtype
#        'keywords': {},  # kwname -> docname, synopsis
    }
    indices = [
#       ChmKwIndex,
    ]

    def clear_doc(self, docname):

        for (typ, name), doc in self.data['objects'].items():
            if doc == docname:
                del self.data['objects'][typ, name]

    def resolve_xref(self, env, fromdocname, builder, typ, target, node,
                     contnode):
        objects = self.data['objects']
        objtypes = self.objtypes_for_role(typ)
        for objtype in objtypes:
            # This is for dealing with :kw: when their :keyword: might have a `:length:` != 4.
            # Corner cases might be nasty..?
            if objtype in (u'keyword', u'command'):
                orig_target = target
                while 1:
                    # look for the target name in the global dictionary, if you dont find it
                    # keep shortening the name you look for until you find it, or fail.
                    # this allows you to have a stylized :kw:, where the first n chars are used
                    # to create the lookup key, and all the characters are printed in the manual
                    if (objtype, target) in objects:
                        # yeah....
                        contnode.children[0] = nodes.Text(CHARmmify(orig_target, len(target)))
                        break
                    else:
                        target = target[:-1]
                    if not target:
                        return
            if (objtype, target) in objects:
                return make_refnode(builder, fromdocname,
                                    objects[objtype, target],
                                    objtype + '-' + target,
                                    contnode, target + ' ' + objtype)


    def get_objects(self):
        for (typ, name), docname in self.data['objects'].iteritems():
            yield name, name, typ, docname, typ + '-' + name, 1



def setup(app):
    app.add_domain(ChmDomain)
