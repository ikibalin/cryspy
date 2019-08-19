"""
define classes to describe crystal 
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
from abc import ABC, abstractmethod


class AbstractCell(ABC):
    """
    Abstract Cell
    """
    @abstractmethod
    def __init__(self, a = 1.0, b = 1.0, c = 1.0, alpha = 90.0, beta = 90.0, 
                 gamma= 90., singony = "Triclinic") -> None:
        pass

    @abstractmethod
    def __repr__(self) -> str:
        pass

    @abstractmethod
    def _refresh(self, a, b, c, alpha, beta, gamma, singony) -> str:
        pass

    @abstractmethod
    def _constr_singony(self) -> None:
        pass

    @abstractmethod
    def set_val(self, a = None, b = None, c = None, alpha = None, 
                   beta = None, gamma= None, singony = None) -> None:
        pass

    @abstractmethod
    def get_val(self, label):
        pass

    @abstractmethod
    def list_vals(self):
        pass

    @abstractmethod
    def _calc_cos_abc(self):
        pass

    @abstractmethod
    def _calc_cos_iabc(self):
        pass

    @abstractmethod
    def _calc_volume(self):
        pass

    @abstractmethod
    def _calc_iucp(self):
        pass

    @abstractmethod
    def _calc_m_b(self):
        pass

    @abstractmethod
    def _calc_m_ib(self):
        pass

    @abstractmethod
    def calc_sthovl(self, h, k, l):
        pass

    @abstractmethod
    def calc_k_loc(self, h, k, l):
        pass

    @abstractmethod
    def calc_m_t(self, h, k, l):
        pass

    @abstractmethod
    def is_variable(self):
        pass

    @abstractmethod
    def get_variables(self):
        pass

    @abstractmethod
    def apply_constraint(self):
        pass

    @abstractmethod
    def calc_hkl(self, space_group, sthovl_min, sthovl_max):
        pass

    @abstractmethod
    def calc_hkl_in_range(self, sthovl_min, sthovl_max):
        pass
