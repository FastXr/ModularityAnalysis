import cobra
from preprocess import ModularityFinder
import json
import os
import operator
import xlsxwriter



class ModularityAnalysis:
    def __init__(self, model_, pathways):
        self.model = cobra.io.load_json_model(model_)
        self.pathway_counts = self.json_converter(pathways)

        self.negative_subsystems = {}
        self.positive_subsystems = {}
        # TODO we may implement module counter
        self.positive_counter = 1
        self.negative_counter = 1


    def report(self, subsystem, name, type_):
        """"
        Input is a dictionary keys are module numbers and values are subsystem dictionaries which
        have systems as a key and how many reactions they have in that module as value and their total number of
        reactions, lastly their percentage coverage
        """
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
        header_list = ['Pathway','Percentage', 'Total', 'in Module', ]
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
    def sorter(dict_):
        # TODO implement function to sort dict according to percentage
        sorted_x = sorted(dict_.items(), key=operator.itemgetter(1), reverse=True)
        return sorted_x


    def analysis(self, dataset_):
        name = self.disease_name(dataset_)
        p = ModularityFinder(self.model)
        dataset = self.json_converter(dataset_)
        positive, negative = p.analysis(dataset_)
        for p_ in positive:  # Top 10
            self.analysis_handler(p_, dataset[0], type_='Positive')
        for n_ in negative:  # Top 10
            self.analysis_handler(n_, dataset[0], type_='Negative')

        #self.sorter(self.positive_subsystems[1])
        self.report(self.positive_subsystems, name, 'Positive')
        self.report(self.negative_subsystems, name, 'Negative')
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

    @staticmethod
    def transport_elimination(dict_):
        black_list = ['Transport', 'Exchange', '_', '']
        to_del = []
        for k in dict_:
            for ban in black_list:
                if k.startswith(ban):
                    to_del.append(k)
        for i in to_del:
            del dict_[i]

        return dict_

    @staticmethod
    def condition_handler(metabolites, modularity):
        for m_ in metabolites:
            if m_.id not in modularity:
                return False
        return True

    @staticmethod
    def json_converter(file):
        with open(file) as f:
            data = json.load(f)
        return data


if __name__ == '__main__':
    m = ModularityAnalysis('./Datasets/recon2.json', './Datasets/pathway_counts.json')
    m.analysis('./Datasets/Enis_bc.json')
    crc = ModularityAnalysis('./Datasets/recon2.json', './Datasets/pathway_counts.json')
    crc.analysis('./Datasets/Enis_CRC.json')
    crohn = ModularityAnalysis('./Datasets/recon2.json', './Datasets/pathway_counts.json')
    crohn.analysis('./Datasets/Enis_crohn.json')