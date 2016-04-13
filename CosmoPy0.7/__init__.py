"""The CosmoPy Module

In this module you will find all classes for theoretical and
statistical cosmological computations.

See the description of individual components for more.

Current revision:
    ID:         $Id: __init__.py 13 2005-09-03 20:30:39Z budavari $
    Date:       $Date: 2005-09-03 10:30:39 -1000 (Sat, 03 Sep 2005) $
    Revision:   $Revision: 13 $
"""

REVISION = '$Revision: 13 $'

def import_all():
    """By default the nested modules are not imported automatically.
    Call this function if you would like to import them all.
    This may be useful for autocompletion in interactive mode."""
    import theory

