"""
Copyright (C) 2008-2011 Volvox Development Team.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Author: Konstantin Lepa <konstantin.lepa@gmail.com>

ANSII Color formatting for output in terminal.
"""

import os

__all__ = ["colored", "cprint"]

VERSION = (1, 1, 0)

ATTRIBUTES = dict(
    list(
        zip(
            ["bold", "dark", "", "underline", "blink", "", "reverse", "concealed"],
            list(range(1, 9)),
        )
    )
)
del ATTRIBUTES[""]

HIGHLIGHTS = dict(
    list(
        zip(
            [
                "on_grey",
                "on_red",
                "on_green",
                "on_yellow",
                "on_blue",
                "on_magenta",
                "on_cyan",
                "on_white",
            ],
            list(range(40, 48)),
        )
    )
)

COLORS = dict(
    list(
        zip(
            [
                "grey",
                "red",
                "green",
                "yellow",
                "blue",
                "magenta",
                "cyan",
                "white",
            ],
            list(range(30, 38)),
        )
    )
)

RESET = "\033[0m"

__ISON = True


def enable(true_false):
    """
    Enable or disable ANSI color formatting.

    Args:
        true_false (bool): True to enable, False to disable.
    """
    global __ISON
    __ISON = true_false


def ison():
    """
    Check if ANSI color formatting is activated.

    Returns:
        bool: True if activated, False otherwise.
    """
    return __ISON


def stream_has_colours(stream):
    """
    Check if the given stream supports colors.

    Args:
        stream: The stream to check (e.g., sys.stdout).

    Returns:
        bool: True if the stream supports colors, False otherwise.
    """
    if not hasattr(stream, "isatty"):
        return False

    if not stream.isatty():
        return False  # auto color only on TTYs
    try:
        import curses

        curses.setupterm()
        return curses.tigetnum("colors") > 2
    except:
        return False  # guess false in case of error


def colored(text, color=None, on_color=None, attrs=None):
    """
    Colorize text for terminal output.

    Args:
        text (str): Text to colorize.
        color (str, optional): Text color.
        on_color (str, optional): Background color.
        attrs (list, optional): List of attributes (e.g., ['bold']).

    Returns:
        str: Colorized text.

    Available text colors:
        red, green, yellow, blue, magenta, cyan, white.

    Available text highlights:
        on_red, on_green, on_yellow, on_blue, on_magenta, on_cyan, on_white.

    Available attributes:
        bold, dark, underline, blink, reverse, concealed.

    Example:
        colored('Hello, World!', 'red', 'on_grey', ['blue', 'blink'])
        colored('Hello, World!', 'green')
    """
    if __ISON and os.getenv("ANSI_COLORS_DISABLED") is None:
        fmt_str = "\033[%dm%s"
        if color is not None:
            text = fmt_str % (COLORS[color], text)

        if on_color is not None:
            text = fmt_str % (HIGHLIGHTS[on_color], text)

        if attrs is not None:
            for attr in attrs:
                text = fmt_str % (ATTRIBUTES[attr], text)

        text += RESET
    return text


def cprint(text, color=None, on_color=None, attrs=None, **kwargs):
    """
    Print colorized text.

    Args:
        text (str): Text to print.
        color (str, optional): Text color.
        on_color (str, optional): Background color.
        attrs (list, optional): List of attributes.
        **kwargs: Arguments passed to the built-in print function.
    """
    try:
        print((colored(text, color, on_color, attrs)), **kwargs)
    except TypeError:
        # flush is not supported by py2.7
        kwargs.pop("flush", None)
        print((colored(text, color, on_color, attrs)), **kwargs)


def colored_map(text, cmap):
    """
    Return colorized text based on a mapping of tokens to colors.

    Args:
        text (str): The text to process.
        cmap (dict): Dictionary mapping tokens to color options.

    Returns:
        str: The colorized text.

    Example:
        colored_key("foo bar", {"bar": "green"})
        colored_key("foo bar", {"bar": {"color": "green", "on_color": "on_red"}})
    """
    if not __ISON:
        return text
    for key, v in cmap.items():
        if isinstance(v, dict):
            text = text.replace(key, colored(key, **v))
        else:
            text = text.replace(key, colored(key, color=v))
    return text


def cprint_map(text, cmap, **kwargs):
    """
    Print colorized text based on a mapping of tokens to colors.

    Args:
        text (str): The text to print.
        cmap (dict): Dictionary mapping tokens to color options.
        **kwargs: Arguments passed to the built-in print function.

    Example:
        cprint_map("Hello world", {"Hello": "red"})
    """
    try:
        print(colored_map(text, cmap), **kwargs)
    except TypeError:
        # flush is not supported by py2.7
        kwargs.pop("flush", None)
        print(colored_map(text, cmap), **kwargs)


def get_terminal_size():
    """
    Return the size of the terminal.

    Returns:
        tuple: (nrow, ncols)

    Based on:
        http://stackoverflow.com/questions/566746/how-to-get-console-window-width-in-python
    """
    try:
        rc = os.popen("stty size", "r").read().split()
        return int(rc[0]), int(rc[1])
    except:
        pass

    env = os.environ

    def ioctl_GWINSZ(fd):
        try:
            import fcntl
            import struct
            import termios

            rc = struct.unpack("hh", fcntl.ioctl(fd, termios.TIOCGWINSZ, "1234"))
            return rc
        except:
            return None

    rc = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)

    if not rc:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            rc = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass

    if not rc:
        rc = (env.get("LINES", 25), env.get("COLUMNS", 80))

    return int(rc[0]), int(rc[1])


if __name__ == "__main__":
    # enable(False)
    print("Current terminal type: %s" % os.getenv("TERM"))
    print("Test basic colors:")
    cprint("Grey color", "grey")
    cprint("Red color", "red")
    cprint("Green color", "green")
    cprint("Yellow color", "yellow")
    cprint("Blue color", "blue")
    cprint("Magenta color", "magenta")
    cprint("Cyan color", "cyan")
    cprint("White color", "white")
    print("-" * 78)

    print("Test highlights:")
    cprint("On grey color", on_color="on_grey")
    cprint("On red color", on_color="on_red")
    cprint("On green color", on_color="on_green")
    cprint("On yellow color", on_color="on_yellow")
    cprint("On blue color", on_color="on_blue")
    cprint("On magenta color", on_color="on_magenta")
    cprint("On cyan color", on_color="on_cyan")
    cprint("On white color", color="grey", on_color="on_white")
    print("-" * 78)

    print("Test attributes:")
    cprint("Bold grey color", "grey", attrs=["bold"])
    cprint("Dark red color", "red", attrs=["dark"])
    cprint("Underline green color", "green", attrs=["underline"])
    cprint("Blink yellow color", "yellow", attrs=["blink"])
    cprint("Reversed blue color", "blue", attrs=["reverse"])
    cprint("Concealed Magenta color", "magenta", attrs=["concealed"])
    cprint("Bold underline reverse cyan color", "cyan", attrs=["bold", "underline", "reverse"])
    cprint("Dark blink concealed white color", "white", attrs=["dark", "blink", "concealed"])
    print("-" * 78)

    print("Test mixing:")
    cprint("Underline red on grey color", "red", "on_grey", ["underline"])
    cprint("Reversed green on red color", "green", "on_red", ["reverse"])

    # Test cprint_keys
    cprint_map("Hello world", {"Hello": "red"})
    cprint_map("Hello world", {"Hello": {"color": "blue", "on_color": "on_red"}})

    # Test terminal size.
    print("terminal size: %s", get_terminal_size())
