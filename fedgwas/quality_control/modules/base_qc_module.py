
from abc import ABC, abstractmethod

class QCModule(ABC):

    def __init__(self) -> None:
        
        self.qc_params = {}
        self.qc_stats = {}
    
    @abstractmethod
    def perform_diagnostic(self, **kwargs):
        pass

    @abstractmethod
    def perform_filtering(self, **kwargs):
        pass

    @abstractmethod
    def perform_diagnostic_incremental(self, **kwargs):
        pass

    @abstractmethod
    def generate_report(self, **kwargs):
        pass