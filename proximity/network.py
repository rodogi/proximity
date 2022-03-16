from proximity import proximity


class Network():
    def __init__(self, graph):
        """A module-agnostic network wrapper
        
        This class creates a network wrapper for graph_tool:graph and
        :networkx:Graph objects.

        Parameters
        ----------
        graph: :graph_tool:`~Graph` or :networkx:`~Graph`

        Raises
        ------
        ValueError if network not an instance of networkx or graph_tool.
        """

        error = ValueError("graph should be an instance of graph_tool.Graph "
                             "or networkx.Graph")
        try:
            module = graph.__module__
        except:
            raise error
        if module == "graph_tool":
            self.module = "gt"
        elif module == "networkx.classes.graph":
            self.module = "nx"
        else:
            raise error
        if graph.__class__.__name__ != 'Graph':
            raise error

        self.Graph = graph

        def get_proximity(self, T, S, n_bins=100, n_iter=1000):
            return proximity(self, T, S, n_bins=n_bins, n_iter=n_iter)