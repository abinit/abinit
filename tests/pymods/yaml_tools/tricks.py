"""
This file solves problems that can only be solved either by a trick or by
a global refactoring such as python 2 compatibility (should be dropped
soon) or problems related to the modified sys.path.
"""


def cstm_isinstance(obj, cls):
    """
    Rely on the true isinstance when possible.

    Only exists because custom sys.path makes issubclass work unexpectedly
    when two different paths are used to access the same class.
    Recursively browses parent classes until it founds a match or reaches
    the top class.

    Args:
        obj: The object to check.
        cls (type or tuple): The class(es) to check against.

    Returns:
        bool: True if obj is an instance of cls, False otherwise.
    """
    if isinstance(obj, cls):
        return True

    if hasattr(obj, "__class__"):
        return cstm_issubclass(obj.__class__, cls)

    return False


def cstm_issubclass(cls_test, cls_ref):
    """
    Rely on the true issubclass when possible.

    Only exists because custom sys.path make issubclass works unexpectedly
    when two different path are used to access the same class.
    Recursively browses parent classes until it founds a match or reaches
    the top class.

    Args:
        cls_test (type): The class to check.
        cls_ref (type or tuple): The class(es) to check against.

    Returns:
        bool: True if cls_test is a subclass of cls_ref, False otherwise.
    """
    if issubclass(cls_test, cls_ref):
        return True

    if isinstance(cls_ref, (tuple, list)):
        for cls in cls_ref:
            if cstm_issubclass(cls_test, cls):
                return True

    elif cls_test.__name__ == cls_ref.__name__:
        return (cls_test.__module__ == "tests." + cls_ref.__module__
                or "tests." + cls_test.__module__ == cls_ref.__module__)

    elif hasattr(cls_test, "__bases__"):
        for cls in cls_test.__bases__:
            if cstm_issubclass(cls, cls_ref):
                return True
    return False
