"""
define classes to describe crystal 
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
from abc import ABC, abstractmethod

class AbstractAtomSite(ABC):
    """
    Abstract Atom Site
    """
    @abstractmethod
    def __repr__(self) -> str:
        pass
