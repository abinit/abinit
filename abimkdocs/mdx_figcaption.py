"""
Definition List Extension for Python-Markdown
=============================================

Added parsing of Figure captions to Python-Markdown.

A simple example:

    ![](http://lorempixel.com/350/150/)
    :   Lorem ipsum dolor sit amet, consectetur adipiscing elit.
        Praesent at consequat magna, faucibus ornare eros. Nam et
        mattis urna. Cras sodales, massa id gravida

Outputs:

    <figure>
        <img alt="" src="http://lorempixel.com/350/150/" />
        <figcaption>
            <p>Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Praesent at consequat magna, faucibus ornare eros. Nam et
            mattis urna. Cras sodales, massa id gravida</p>
        </figcaption>
    </figure>


Copyright 2013 - [Helder Correia](http://heldercorreia.com)

"""

from __future__ import unicode_literals
from markdown import Extension
from markdown.inlinepatterns import IMAGE_LINK_RE, IMAGE_REFERENCE_RE
from markdown.blockprocessors import BlockProcessor
from markdown.util import etree
import re

import logging
logger = logging.getLogger('MARKDOWN')

FIGURES = [IMAGE_LINK_RE, IMAGE_REFERENCE_RE]


class FigcaptionProcessor(BlockProcessor):
    """ Process figure captions."""

    RE = re.compile(r'(^|\n)[ ]{0,3}:[ ]{1,3}(?P<caption>.*?)(\n|$)')
    FIGURES_RE = re.compile('|'.join(f for f in FIGURES))
    NO_INDENT_RE = re.compile(r'^[ ]{0,3}[^ :]')

    def test(self, parent, block):
        return bool(self.RE.search(block))

    def run(self, parent, blocks):
        # pop the entire block as a single string
        raw_block = blocks.pop(0)

        # Get list of figure elements before the colon (:)
        m = self.RE.search(raw_block)

        # Get elements
        elements = raw_block[:m.start()]
        test_elements = self.FIGURES_RE.search(elements)

        if not test_elements:
            # This is not a figure item.
            blocks.insert(0, raw_block)
            return False

        # Get caption
        block = raw_block[m.end():]
        no_indent = self.NO_INDENT_RE.match(block)

        if no_indent:
            caption, theRest = (block, None)
        else:
            caption, theRest = self.detab(block)
        if caption:
            caption = '%s\n%s' % (m.group('caption'), caption)
        else:
            caption = m.group('caption')

        # Create figure
        figure = etree.SubElement(parent, 'figure')
        figure.text = elements

        # Add definition
        self.parser.state.set('fig')
        figcaption = etree.SubElement(figure, 'figcaption')
        self.parser.parseBlocks(figcaption, [caption])
        self.parser.state.reset()

        if theRest:
            blocks.insert(0, theRest)


class FigcaptionExtension(Extension):
    """ Add definition lists to Markdown. """

    def extendMarkdown(self, md, md_globals):
        """ Add an instance of FigcaptionProcessor to BlockParser. """

        # def_list = 'def_list' in md.registeredExtensions
        md.parser.blockprocessors.add('figcaption',
                                      FigcaptionProcessor(md.parser),
                                      '<ulist')


def makeExtension(**kwargs):
    return FigcaptionExtension(**kwargs)
