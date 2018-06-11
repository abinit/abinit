# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

from mkdocs.plugins import BasePlugin


class AbiMkdocsPlugin(BasePlugin):

    #def on_config(self, config, **kwargs):
    #    config['theme'].static_templates.add('my_template.html')
    #    return config

    def on_page_markdown(self, markdown, **kwargs):
        """
        The page_markdown event is called after the page's markdown is loaded from file
        and can be used to alter the Markdown source text.
        The meta-data has been stripped off and is available as page.meta at this point.

        Args:
            markdown: Markdown source text of page as string
            page: mkdocs.nav.Page instance
            config: global configuration object
            site_navigation: global navigation object

        Returns:
            Markdown source text of page as string
        """
        page = kwargs["page"]
        #print("page", page, "\ninput_path", page.input_path, "\npage.meta", page.meta)
        #if "authors" in page.meta:
        #    print("authors:", type(page.meta["authors"]), page.meta["authors"])

        #config = kwargs["config"]
        #for key, value in config["abimkdocs_links"].items():
        #    markdown = markdown.replace(key, value)


        # (Re)Add header with metadata and rpath
        mymeta = page.meta.copy()
        mymeta["rpath"] = "/" + page.input_path

        header = "---\n" + "\n".join("%s: %s" % item for item in mymeta.items()) + "\n---\n\n"
        #print(header)

        return header + markdown
