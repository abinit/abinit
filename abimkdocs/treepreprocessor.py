from __future__ import print_function, division, unicode_literals, absolute_import

from markdown import Extension
from markdown.treeprocessors import Treeprocessor


class AbinitTreeprocessor(Treeprocessor):

    def __init__(self, site_navigation, strict):
        self.site_navigation = site_navigation
        self.strict = strict

    def set_link_class(self, element):
        for child in element:
            if child.tag == "a":
                child.set("class", "myclass") #set the class attribute
            set_link_class(child) # run recursively on children

    def run(self, root):
        """Update urls on anchors and images to make them relative
        Iterates through the full document tree looking for specific
        tags and then makes them relative based on the site navigation
        """

        for element in root.iter():

            if element.tag == 'a':
                key = 'href'
            elif element.tag == 'img':
                key = 'src'
            else:
                continue

            url = element.get(key)
            new_url = path_to_url(url, self.site_navigation, self.strict)
            element.set(key, new_url)

        return root


class AbinitTreePreprocessorExtension(Extension):
    """
    The Extension class is what we pass to markdown, it then
    registers the Treeprocessor.
    """

    def __init__(self, site_navigation, strict):
        self.site_navigation = site_navigation
        self.strict = strict

    def extendMarkdown(self, md, md_globals):
        relpath = AbinitTreeprocessor(self.site_navigation, self.strict)
        md.treeprocessors.add("relpath", relpath, "_end")


def makeExtension(*args, **kwargs):
    return AbinitTreeprocessorExtension(*args, **kwargs)
