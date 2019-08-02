import cobra
from preprocess import ModularityFinder
import json
import operator
# import itertools
import xlsxwriter
import pickle
import datetime
import graphviz as gv
from scipy.stats import hypergeom

# TODO Transport Elimination can be done elsewhere


class ModularityAnalysis:
    def __init__(self, model_, pathways, currency_metabolites):
        self.model = cobra.io.load_json_model(model_)
        self.pathway_counts = self.json_converter(pathways)
        self.currency_metabolites = self.json_converter(currency_metabolites)

        self.positive_reactions = {}
        self.negative_reactions = {}

        self.positive_modules = {}
        self.negative_modules = {}
        self.positive_connections = {}
        self.negative_connections = {}



        self.negative_subsystems = {}
        self.positive_subsystems = {}
        self.positive_counter = 1
        self.negative_counter = 1

    def analysis_new(self, dataset_, scale=None):
        name = self.disease_name(dataset_)
        if scale != None:
            name += '_scaled{}_{}'.format(scale[0], scale[1])
        p = ModularityFinder(self.model)
        dataset = self.json_converter(dataset_)
        positive, negative = p.analysis(dataset_, scale=scale)

        print('New analysis started')
        print(datetime.datetime.now())

        positive_depth, negative_depth = p.network_sizes()
        positive_generator = self.generator_iterator(positive, positive_depth)
        positive_updated = tuple(c for c in next(positive_generator))

        for p_ in positive_updated:  # Top 10
            self.analysis_handler(p_, dataset[0], type_='Positive')

        print('Report Part')
        print(datetime.datetime.now())
        self.connection_analysis(self.positive_modules, type_='Positive')
        self.enrichment_analysis(self.positive_reactions, name + '_depth{}'.format(positive_depth), 'Positive')



        print('Negative Part Stars')
        print(datetime.datetime.now())

        negative_generator = self.generator_iterator(negative, negative_depth)

        negative_updated = tuple(c for c in next(negative_generator))

        for p_ in negative_updated:  # Top 10
            self.analysis_handler(p_, dataset[0], type_='Negative')

        self.connection_analysis(self.negative_modules, type_='Negative')
        self.enrichment_analysis(self.negative_reactions, name + '_depth{}'.format(negative_depth), 'Negative')

        print(datetime.datetime.now())
        print('END OF ALL')

    def enrichment_analysis(self, module_, name, type_):
        """
        :param module_: Dictionary keys are modules values are reactions in that module
        :return:
        """
        print('ENRICHMENT ANALYSIS')
        print(module_)  #  This print statement can be used to implement functions in jupyter

        total_dict = {}
        module_dict = {}
        result_dict = {}

        for m_ in module_:
            module_dict.setdefault(m_, {})
            for r_ in module_[m_]:
                s = self.model.reactions.get_by_id(r_).subsystem

                total_dict.setdefault(s, set())
                total_dict[s].add(r_)

                module_dict[m_].setdefault(s, set())
                module_dict[m_][s].add(r_)


        N = sum(map(len, total_dict.values()))


        for m_ in module_dict:
            for reaction in module_dict[m_]:
                k = len(module_dict[m_][reaction])
                n = len(module_dict[m_])
                K = len(total_dict[reaction])

                # r = k / (n * K / N)
                p = hypergeom.pmf(k=k, M=N, n=K, N=n)
                if p < 0.05: # Enrichment Analysis Threshold  # TODO
                    result_dict.setdefault(m_, {})
                    result_dict[m_].setdefault(reaction, [])
                    result_dict[m_][reaction] = [len(module_dict[m_][reaction]), len(total_dict[reaction]), p]

        self.graph_draw(result_dict, name=name, type_=type_)
        self.report(result_dict, name=name, type_=type_)

    def analysis_handler(self, module, dataset, type_='Positive'):
        # That c variable is being to used to determine which module we are working on can be updated.
        if type_ == 'Positive':
            subsystem = self.positive_subsystems
            reactions = self.positive_reactions
            modules_ = self.positive_modules
            c = self.positive_counter
            self.positive_counter += 1
        elif type_ == 'Negative':
            subsystem = self.negative_subsystems
            reactions = self.negative_reactions
            modules_ = self.negative_modules
            c = self.negative_counter
            self.negative_counter += 1

        module = self.module_updater(module)  # 03.07 To get rid of currency metabolites in a module
        modules_[c] = module
        possible = self.reaction_needs(module, dataset)
        reactions[c] = possible
        to_save = self.subsystems_calculator(possible)
        subsystem[c] = to_save

    def subsystems_calculator(self, modularity):
        # This function determines which subsystems are taking part in certain module so that we can try to determine whether they are enrich or not.
        """
        :param modularity: Reactions in a module
        :return: Will return an enriched analysed thing
        """


        dict_to_return = {}
        for i in modularity:
            s = self.model.reactions.get_by_id(i).subsystem
            dict_to_return.setdefault(s, [self.pathway_counts[s], 0])
            dict_to_return[s][1] += 1
        for system in dict_to_return:
            t = dict_to_return[system][0]  # Total
            m = dict_to_return[system][1]  # in Module
            percent = (m*100.0)/t  # Percentage
            dict_to_return[system] = [percent, t, m]

        return dict_to_return

    def reaction_needs(self, modularity_analysis, reactions):
        possible_reactions = {}
        for metabolite in modularity_analysis:
            for r in self.model.metabolites.get_by_id(metabolite).reactions:
                if self.condition_handler(self.model.reactions.get_by_id(r.id).metabolites, modularity_analysis):
                    possible_reactions.setdefault(r.id, set()) # To debug
                    possible_reactions[r.id].add(reactions[r.id])
        return possible_reactions

    def report(self, subsystem, name, type_):
        """"
        Input is a dictionary keys are module numbers and values are subsystem dictionaries which
        have systems as a key and how many reactions they have in that module as value and their total number of
        reactions, lastly their percentage coverage
        """
        subsystem = self.transport_elimination(subsystem)
        """
        Pickling eliminated 12.06
        """
        # with open('./Results/{}__{}.pickle'.format(name, type_), 'wb') as outfile:
        #     pickle.dump(subsystem, outfile)
        print(subsystem)
        workbook = xlsxwriter.Workbook('./Results/{}_{}.xlsx'.format(name, type_))
        if type_ == 'Positive':
            sheet = workbook.add_worksheet('Positive')
        elif type_ == 'Negative':
            sheet = workbook.add_worksheet('Negative')
        header_list = ['Pathway', 'In Module', 'Total', 'P Value', ]
        row = 0
        column = 1
        for i in header_list:
            sheet.write(row, column, i)
            column += 1
        row = 1
        for n_ in subsystem:
            column = 0
            sheet.write(row, column, n_)
            column += 1
            row += 1
            iter_ = self.sorter(subsystem[n_])
            for i in iter_:
                column = 1
                sheet.write(row, column, i[0])
                # Filter to get ride of below 10 percent
                # if i[1][0] < 10:  # Filter is removed 09/03/2019
                #     continue
                for v_ in i[1]:
                    column += 1
                    sheet.write(row, column, v_)
                row += 1
        workbook.close()

    def module_updater(self, module):
        return [metabolite for metabolite in module if metabolite not in self.currency_metabolites]


    def graph_draw(self, dict_, name, type_):
        """
        :param dict_: dictionary of clusters
        :param name: naming for the file
        :return:
        This function only creates proper pictures

        Added at 12.06 May need to be changed after Enrichment Analysis
        """
        existing_clusters = set()
        g = gv.Digraph('G',
                       filename='./Graphs/_{name}.gv'.format(name=name),
                       engine='fdp')
        for cluster in dict_:
            if len(dict_[cluster]) < 2:
                continue

            with g.subgraph(name='cluster_{}'.format(cluster)) as c:
                existing_clusters.add(cluster)
                c.attr(label='cluster_{}'.format(cluster))
                for node_ in dict_[cluster]:
                    c.node(node_)
        if type_ == 'Positive':
            connection_ = self.positive_connections
        elif type_ == 'Negative':
            connection_ = self.negative_connections

        for source in connection_:
            for destination in connection_[source]:
                if source in existing_clusters and destination in existing_clusters:
                    g.edge('cluster_{}'.format(source), 'cluster_{}'.format(destination),
                           len='20.0', width='2.0', label="{}".format(connection_[source][destination]))


        g.save()


    def connection_analysis(self, module_, type_):
        connection_analysis = {}
        for clust_ in module_:
            for metabolite in module_[clust_]:
                for r_ in self.model.metabolites.get_by_id(metabolite).reactions:
                    input_metabolites = set()
                    output_metabolites = set()
                    meta = self.model.reactions.get_by_id(r_.id).metabolites
                    for m_ in meta:
                        if meta[m_] == 1.0:
                            output_metabolites.add(m_.id)
                        if meta[m_] == -1.0:
                            input_metabolites.add(m_.id)
                    if input_metabolites.issubset(set(module_[clust_])):
                        for other_clust in module_:
                            if clust_ == other_clust or output_metabolites == set():
                                continue
                            if output_metabolites.issubset(set(module_[other_clust])):
                                connection_analysis.setdefault(clust_, {})
                                connection_analysis[clust_].setdefault(other_clust, 0)
                                connection_analysis[clust_][other_clust] += 1
        if type_=='Positive':
            self.positive_connections = connection_analysis
        elif type_=='Negative':
            self.negative_connections = connection_analysis

    @staticmethod
    def transport_elimination(dict_):
        copy_dict = dict_
        black_list = ['Transport', 'Exchange', '_']
        to_del = []
        for m in dict_:
            for k in dict_[m]:
                for ban in black_list:
                    if k.startswith(ban):
                        to_del.append((m, k))
        for p in to_del:
            del copy_dict[p[0]][p[1]]

        return copy_dict

    @staticmethod
    def condition_handler(metabolites, modularity):
        for m_ in metabolites:
            if m_.id not in modularity:
                return False
        return True

    @staticmethod
    def sorter(dict_):
        return sorted(dict_.items(), key=operator.itemgetter(1), reverse=True)

    @staticmethod
    def json_converter(file):
        with open(file) as f:
            data = json.load(f)
        return data

    @staticmethod
    def generator_iterator(generator, n):
        for i in range(n-1):
            next(generator)
        return generator

    @staticmethod
    def disease_name(str_):
        l_ = str_.split('/')
        name_with_ext = l_[-1]
        n = name_with_ext.split('.')
        name = n[0]
        return name
if __name__ == '__main__':
    m = ModularityAnalysis('./Datasets/recon2.json',
                           './Datasets/pathway_counts.json',
                           './Datasets/currency_metabolites.json')

    m.analysis_new('./Datasets/Enis_bc.json', scale=(-1, 1))

    # m.report(r, 'EnisSon', type_='Positive')




    """
    scale_list = [(-1, 1), (-10, 10), (-100, 100)]
    dataset_list = ['./Datasets/Enis_bc.json',
                    './Datasets/Enis_CRC.json',
                    './Datasets/Enis_crohn.json']

    for d_ in dataset_list:
        for scale_ in scale_list:
            m = ModularityAnalysis('./Datasets/recon2.json', './Datasets/pathway_counts.json')
            m.analysis_new(d_, scale=scale_)

    """


