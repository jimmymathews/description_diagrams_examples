
from .diagram_mutation import DescriptionDiagramMutator


class DiagramAnnotator(DescriptionDiagramMutator):
    def mutate(self, diagram):
        self.annotate(diagram)

    def annotate(self, diagram):
        for s1, s2 in diagram.graph.edges:
            l = diagram.graph.edges[s1, s2]['supports feature inference']
            diagram.graph.edges[s1, s2]['supports feature inference'] = 0
            diagram.graph.edges[s1, s2]['number of feature inferences supported'] = len(l)
            diagram.graph.edges[s1, s2]['object set size'] = len(set([obj for obj, feature in l]))
            diagram.graph.edges[s1, s2]['feature set size'] = len(set([feature for obj, feature in l]))
        for s in diagram.graph.nodes:
            diagram.graph.nodes[s]['number of features'] = len(s)
            diagram.graph.nodes[s]['binary representation'] = diagram.facet(s).get_binary_vector_string()
