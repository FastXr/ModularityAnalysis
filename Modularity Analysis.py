import cobra
from preprocess import ModularityFinder
import json
import operator
import itertools
import xlsxwriter

import datetime

class ModularityAnalysis:
    def __init__(self, model_, pathways):
        self.model = cobra.io.load_json_model(model_)
        self.pathway_counts = self.json_converter(pathways)

        self.negative_subsystems = {}
        self.positive_subsystems = {}
        self.positive_counter = 1
        self.negative_counter = 1

    def analysis_new(self, dataset_, scale=None, depth=1):
        name = self.disease_name(dataset_)
        # TODO add depth to name if you exceed depth
        if scale != None:
            name += '_scaled{}_{}'.format(scale[0], scale[1])


        p = ModularityFinder(self.model)
        dataset = self.json_converter(dataset_)
        positive, negative = p.analysis(dataset_, scale=scale)

        print('New anaylsis started')
        print(datetime.datetime.now())


        # Changed by me at 25/02 20.30 put it back This is working BTW
        positive_generator = self.generator_iterator(positive, depth)



        #limited = itertools.takewhile(lambda c: len(c) <= depth, positive)
        #positive_generator = self.generator_iterator(limited, depth-2)



        print('Generator Ended')

        #positive_generator = self.generator_iterator(limited, n=int(depth/2))

        positive_updated = tuple(c for c in next(positive_generator))
        #positive_updated = tuple(c for c in positive_generator)
        #print(positive_updated)
        for p_ in positive_updated:  # Top 10
            self.analysis_handler(p_, dataset[0], type_='Positive')
        print('Report Part')
        print(datetime.datetime.now())
        self.report(self.positive_subsystems, name + '_depth{}'.format(depth), 'Positive')
        print(datetime.datetime.now())
        print('Ends')
        # TODO ADD NEGATIVE PART IF WORKS

        print('Negative Part Stars')
        print(datetime.datetime.now())

        negative_generator = self.generator_iterator(negative, depth)

        negative_updated = tuple(c for c in next(negative_generator))

        for p_ in negative_updated:  # Top 10
            self.analysis_handler(p_, dataset[0], type_='Negative')
        self.report(self.negative_subsystems, name + '_depth{}'.format(depth), 'Negative')

        print(datetime.datetime.now())
        print('END OF ALLLL')












    def analysis_old(self, dataset_, scale=None, depth=1):
        name = self.disease_name(dataset_)
        # TODO add depth to name if you exceed depth
        if scale != None:
            name += '_scaled{}_{}'.format(scale[0], scale[1])

        p = ModularityFinder(self.model)
        dataset = self.json_converter(dataset_)
        positive, negative = p.analysis(dataset_, scale=scale)

        # TODO Update This part since our dataset changed
        k = depth
        print('Last Part Started')

        # positive_updated = tuple(c for c in next(positive))


        for communities in itertools.islice(positive, k):
            positive_updated = tuple(c for c in communities)
        # for communities in itertools.islice(negative, k):
        #     negative_updated = tuple(sorted(c) for c in communities)
        print('Part middle')
        for p_ in positive_updated:  # Top 10
            self.analysis_handler(p_, dataset[0], type_='Positive')
        print('Last part Ended')


        self.report(self.positive_subsystems, name, 'Positive')
        # self.report(self.negative_subsystems, name, 'Negative')

    def analysis(self, dataset_, scale=None, depth=1):
        name = self.disease_name(dataset_)
        # TODO add depth to name if you exceed depth
        if scale != None:
            name += '_scaled{}_{}'.format(scale[0],scale[1])

        p = ModularityFinder(self.model)
        dataset = self.json_converter(dataset_)
        positive, negative = p.analysis(dataset_, scale=scale)

        # TODO Update This part since our dataset changed
        k = depth
        print('Come HERE')
        depth_positive = 1
        for communities in itertools.islice(positive, k):
            print('ENISSSS')
            print('Positive{}'.format(depth_positive))
            self.negative_subsystems = {}
            self.positive_subsystems = {}
            self.positive_counter = 1
            self.negative_counter = 1
            # TODO we may use for loop to iterate to curtain point consider it. !!!!!
            positive_updated = tuple(c for c in communities)
            for p_ in positive_updated:  # Top 10
                self.analysis_handler(p_, dataset[0], type_='Positive')
            self.report(self.positive_subsystems, name+'_depth{}'.format(depth_positive), 'Positive')
            print(datetime.datetime.now())
            depth_positive += 1

        depth_negative = 1
        for communities in itertools.islice(negative, k):
            """
            This part is added by stealing from __init__
            We may change location 
            """
            print('Negative'.format(depth_negative))
            self.negative_subsystems = {}
            self.positive_subsystems = {}
            self.positive_counter = 1
            self.negative_counter = 1
            negative_updated = tuple(c for c in communities)
            for n_ in negative_updated:  # Top 10
                self.analysis_handler(n_, dataset[0], type_='Negative')
            self.report(self.negative_subsystems, name+'_depth{}'.format(depth_negative), 'Negative')
            print(datetime.datetime.now())
            depth_negative += 1



        # for p_ in positive:  # Top 10
        #     self.analysis_handler(p_, dataset[0], type_='Positive')
        # for n_ in negative:  # Top 10
        #     self.analysis_handler(n_, dataset[0], type_='Negative')

        #self.sorter(self.positive_subsystems[1])



        #print(self.positive_subsystems)
        #print(self.negative_subsystems)
        #print(self.positive_counter)
        #print(self.negative_counter)


        # print(sorted(self.positive_subsystems.items(), key=operator.itemgetter(1), reverse=True))
        # print(sorted(self.negative_subsystems.items(), key=operator.itemgetter(1), reverse=True))

    @staticmethod
    def disease_name(str_):
        l_ = str_.split('/')
        name_with_ext = l_[-1]
        n = name_with_ext.split('.')
        name = n[0]
        return name

    def analysis_handler(self, module, dataset, type_='Positive'):
        if type_ == 'Positive':
            subsystem = self.positive_subsystems
            c = self.positive_counter
            self.positive_counter += 1
        elif type_ == 'Negative':
            subsystem = self.negative_subsystems
            c = self.negative_counter
            self.negative_counter += 1

        possible = self.reaction_needs(module, dataset)
        to_save = self.subsystems_calculator(possible)
        subsystem[c] = to_save

    def subsystems_calculator(self, modularity):
        dict_to_return = {}
        for i in modularity:
            s = self.model.reactions.get_by_id(i).subsystem
            dict_to_return.setdefault(s, [self.pathway_counts[s], 0])
            dict_to_return[s][1] += 1
        for system in dict_to_return:
            t = dict_to_return[system][0]  # Total
            m = dict_to_return[system][1]  # in Module
            percent = (m*100)/t  # Percentage
            dict_to_return[system] = [percent, t, m]
        return dict_to_return

    def reaction_needs(self, modularity_analysis, reactions):
        possible_reactions = {}
        for metabolite in modularity_analysis:
            for r in self.model.metabolites.get_by_id(metabolite).reactions:
                # if r.id in reactions:
                if self.condition_handler(self.model.reactions.get_by_id(r.id).metabolites, modularity_analysis):
                    possible_reactions.setdefault(r.id, [])   # To debug
                    if reactions[r.id] not in possible_reactions[r.id]:
                        possible_reactions[r.id].append(reactions[r.id])
        return possible_reactions
    def report(self, subsystem, name, type_):
        """"
        Input is a dictionary keys are module numbers and values are subsystem dictionaries which
        have systems as a key and how many reactions they have in that module as value and their total number of
        reactions, lastly their percentage coverage
        """

        subsystem = self.transport_elimination(subsystem)
        print(subsystem)
        # TODO Check if file exist then act accordingly
        workbook = xlsxwriter.Workbook('./Results/{}_{}.xlsx'.format(name, type_))
        # IF necessary
        # TODO implement if condition and type_
        # TODO implement different sheets
        if type_ == 'Positive':
            sheet = workbook.add_worksheet('Positive')
        elif type_ == 'Negative':
            sheet = workbook.add_worksheet('Negative')
        header_list = ['Pathway', 'Percentage', 'Total', 'in Module', ]
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
            iter = self.sorter(subsystem[n_])
            for i in iter:
                column = 1
                sheet.write(row, column, i[0])
                for v_ in i[1]:
                    column += 1
                    sheet.write(row, column, v_)
                row += 1

            # for s_ in subsystem[n_]:
            #     column = 1
            #     sheet.write(row, column, s_)
            #     # For loop start
            #     for v_ in subsystem[n_][s_]:
            #         column += 1
            #         sheet.write(row, column, v_)
            #     # For loop end
            #     row += 1
        workbook.close()
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
        # TODO implement function to sort dict according to percentage
        sorted_x = sorted(dict_.items(), key=operator.itemgetter(1), reverse=True)
        return sorted_x
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

if __name__ == '__main__':
    # m = ModularityAnalysis('./Datasets/recon2.json', './Datasets/pathway_counts.json')
    # m.analysis('./Datasets/Enis_bc.json', scale=(-1, 1))
    #
    # m = ModularityAnalysis('./Datasets/recon2.json', './Datasets/pathway_counts.json')
    # m.analysis('./Datasets/Enis_bc.json', scale=(-10, 10))
    #
    # m = ModularityAnalysis('./Datasets/recon2.json', './Datasets/pathway_counts.json')
    # m.analysis('./Datasets/Enis_bc.json', scale=(-100, 100))

    m = ModularityAnalysis('./Datasets/recon2.json', './Datasets/pathway_counts.json')
    #m.analysis_old('./Datasets/Enis_bc.json', depth=1)
    #m.analysis_new('./Datasets/Enis_bc.json', depth=7)
    #m.analysis('./Datasets/Enis_bc.json', depth=5)

    #m = ModularityAnalysis('./Datasets/recon2.json', './Datasets/pathway_counts.json')
    #m.analysis_new('./Datasets/Enis_bc.json', depth=25, scale=(-1, 1))

    m = ModularityAnalysis('./Datasets/recon2.json', './Datasets/pathway_counts.json')
    m.analysis_new('./Datasets/Enis_bc.json', depth=35, scale=(-1, 1))




    # crc = ModularityAnalysis('./Datasets/recon2.json', './Datasets/pathway_counts.json')
    # crc.analysis('./Datasets/Enis_CRC.json')
    # crohn = ModularityAnalysis('./Datasets/recon2.json', './Datasets/pathway_counts.json')
    # crohn.analysis('./Datasets/Enis_crohn.json')