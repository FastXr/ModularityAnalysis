import networkx as nx
import cobra
import json
from math import sqrt


class ModularityFinder:
    def __init__(self, model_):
        self.positive_network = nx.MultiDiGraph()
        self.negative_network = nx.MultiDiGraph()
        try:
            self.model = cobra.io.load_json_model(model_)
        except:
            self.model = model_

    def analysis(self, reactions, scale=None):
        self.case_network_builder(reactions, scale=scale)
        positive_analysis = nx.algorithms.community.centrality.girvan_newman(self.positive_network)
        negative_analysis = nx.algorithms.community.centrality.girvan_newman(self.negative_network)
        return positive_analysis, negative_analysis

    def network_sizes(self):
        return int(sqrt(len(self.positive_network))), int(sqrt(len(self.negative_network)))

    def case_network_builder(self, reactions, scale=None):
        # Scale is a tuple that eliminates values between those values
        reactions_ = self.json_converter(reactions)
        reactions__ = self.average_score_calculator(reactions_)
        for reaction in reactions__:
            v = reactions__[reaction]
            if scale != None:
                if scale[0] < v < scale[1]:
                    continue
            if v > 0:
                self.add_edge(reaction, v, self.positive_network)
            elif v <= 0:
                self.add_edge(reaction, v, self.negative_network)

    def add_edge(self, reaction, value, network):
        if value < 0:
            value = -value
        metabolites_ = self.model.reactions.get_by_id(reaction).metabolites.items()
        input_metabolites = []
        output_metabolites = []
        for metabolite_pair in metabolites_:
            metabolite = metabolite_pair[0].id
            network.add_node(metabolite)
            if metabolite_pair[1] == -1.0:
                input_metabolites.append(metabolite)
            else:
                output_metabolites.append(metabolite)
        for input_ in input_metabolites:
            for out_ in output_metabolites:
                network.add_edge(input_, out_, weight=value)

    @staticmethod
    def average_score_calculator(reaction_values):
        """
        Reaction values is a list of dictionaries
        Dictionaries are different people's results.
        Preprocess part
        """
        average_scoring_dict = {}
        for person_dict in reaction_values:
            for reaction in person_dict:
                average_scoring_dict.setdefault(reaction, 0)
                average_scoring_dict[reaction] += person_dict[reaction]
        for r in average_scoring_dict:
            prev_value = average_scoring_dict[r]
            prev_value = prev_value / len(reaction_values)
            average_scoring_dict[r] = prev_value
        return average_scoring_dict

    def sample_network_builder(self):
        """
        :return:
        This function was used to discover needs for this project as building a decoy network
        """
        for i in range(len(self.model.reactions)):
            reaction = self.model.reactions.__getitem__(i)
            metabolites = reaction.metabolites.items()
            input_metabolites = []
            output_metabolites = []
            for metabolite_pair in metabolites:
                metabolite = metabolite_pair[0].id
                self.positive_network.add_node(metabolite)
                if metabolite_pair[1] == -1.0:
                    input_metabolites.append(metabolite)
                else:
                    output_metabolites.append(metabolite)
            for input_ in input_metabolites:
                for out_ in output_metabolites:
                    self.positive_network.add_edge(input_, out_, weight=1.0)

    @staticmethod
    def json_converter(file):
        with open(file) as f:
            data = json.load(f)
        return data


if __name__ == '__main__':
    m = ModularityFinder('./Datasets/recon2.json')
    m.analysis('./Datasets/Enis_bc.json')