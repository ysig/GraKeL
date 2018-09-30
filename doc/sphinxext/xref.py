# -*- coding: utf-8 -*-
# Taken from: https://github.com/michaeljones/sphinx-xref/blob/master/xref.py

from docutils import nodes

from sphinx.util import caption_ref_re

def xref( typ, rawtext, text, lineno, inliner, options={}, content=[] ):

    title = target = text
    titleistarget = True
    # look if explicit title and target are given with `foo <bar>` syntax
    brace = text.find('<')
    if brace != -1:
        titleistarget = False
        m = caption_ref_re.match(text)
        if m:
            target = m.group(2)
            title = m.group(1)
        else:
            # fallback: everything after '<' is the target
            target = text[brace+1:]
            title = text[:brace]

    link = xref.links[target]

    if brace != -1:
        pnode = nodes.reference(target, title, refuri=link[1])
    else:
        pnode = nodes.reference(target, link[0], refuri=link[1])

    return [pnode], []

def get_refs(app):

    xref.links = app.config.xref_links

def setup(app):

    app.add_config_value('xref_links', {}, True)
    app.add_role('xref', xref)
    app.connect("builder-inited", get_refs)

