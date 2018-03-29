# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

from mkdocs.plugins import BasePlugin


class AbiMkdocsPlugin(BasePlugin):

    #def on_config(self, config, **kwargs):
    #    config['theme'].static_templates.add('my_template.html')
    #    return config

    def on_page_markdown(self, markdown, **kwargs):
        """
        markdown: Markdown source text of page as string
        page: mkdocs.nav.Page instance
        config: global configuration object
        site_navigation: global navigation object

        Returns:
            Markdown source text of page as string
        """
        page = kwargs["page"]
        #print("page", page)
        #print("input_path", page.input_path)
        #print("page.meta", page.meta)
        #if "authors" in page.meta:
        #    print("authors:", type(page.meta["authors"]), page.meta["authors"])

        mymeta = page.meta.copy()
        mymeta["rpath"] = "/" + page.input_path

        header = "---\n" + "\n".join("%s: %s" % item for item in mymeta.items())  + "\n---\n\n"
        #print(header)
        s =  header + markdown
        return s
