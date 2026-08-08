from abc import ABC, abstractmethod

class Loader(ABC):
    @abstractmethod
    def load_gene(self, gene_id: str):
        pass

    @abstractmethod
    def gene_list(self):
        pass

    @abstractmethod
    def cache_id(self):
        pass

class FileBasedLoader(Loader):
    def __init__(self):
        self.files_loaded = False

    @abstractmethod
    def load_files(self):
        # cache_id() may be executed before files are loaded
        pass

    def load_gene(self):
        if not self.files_loaded:
            self.load_files()
            self.files_loaded = True

    def gene_list(self):
        if not self.files_loaded:
            self.load_files()
            self.files_loaded = True

    def delete(self):
        pass

    def __del__(self):
        self.delete()