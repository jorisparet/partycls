"""
Useful functions for an internal use in the main code.
"""

import re
import enum
from partycls.particle import aliases


def tipify(s):
    """
    Convert a string ``s`` into the best matching type,
    *i.e.* an instance of ``int``, ``float``, ``str`` 
    or a list of those types.

    Parameters
    ----------
    s : str
        String to convert.

    Returns
    -------
    s : int, float, str or list
        Best-matching type for the input string ``s``.
    """
    if '_' in s:
        return s
    elif ',' in s:
        s = s.split(sep=',')
        s = [tipify(s_i) for s_i in s]
        return s
    else:
        try:
            return int(s)
        except ValueError:
            try:
                return float(s)
            except ValueError:
                return s


def standardize_condition(condition):
    """
    Check that the condition is correctly formated 
    (*i.e* in the form ``"<attr> _operator_ <val>"``).

    Parameters
    ----------
    condition : str
        Condition to standardize.

    Raises
    ------
    ValueError
        If ``condition`` is not valid or if the 
        ``"<attr>"`` keyword is not recognized.

    Returns
    -------
    condition : str
        The standardized condition.
    """
    regexp = re.search('(\w+\.?\w*\[?\d?\]?)\s?(<|<=|\!=|==|>=|>)\s?([\'|\"]?[\w\(\),]+[\'|\"]?)', condition)
    if regexp:
        attr = regexp.group(1)
        operator = regexp.group(2)
        value = regexp.group(3)
        # if attribute is an alias, replace it with full attribute
        if not attr.startswith('particle'):
            if attr in aliases:
                attr = aliases[attr]
            else:
                raise ValueError('attribute "{}" is not recognized'.format(attr))
        # reformat condition
        condition = '{} {} {}'.format(attr, operator, value)
        return condition
    else:
        raise ValueError('"{}" is not a valid condition'.format(condition))

# Valid nearest neighbors methods
_nearest_neighbors_methods_ = ['auto', 'fixed', 'sann', 'voronoi']