"""
define classes to describe crystal 
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
from abc import ABC, abstractmethod

class AbstractFract(ABC):
    """
    Abstract Fract
    """
    @abstractmethod
    def __repr__(self) -> str:
        pass
