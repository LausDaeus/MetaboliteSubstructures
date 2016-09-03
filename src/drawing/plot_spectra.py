import numpy as np
import matplotlib.pyplot as plt
import src.utils.fileIO as io
import sys

term_dict = {"energy2":"high","energy1":"medium","energy0":"low"}#CFM has a shoddy relabeling practice.

def plot_spectra(spectra_dict):
    for spectra in spectra_dict:
        mass_intensities = spectra_dict[spectra]
        if mass_intensities:
            mass_intensities = np.array(mass_intensities)
            plt.bar(mass_intensities[:,0], mass_intensities[:,1])
            plt.ylabel("intensity")
            plt.xlabel("mass/charge")
            plt.title(spectra)
            plt.show()


def plot_cfm_data(cfm_dict, spectra_dict, output_path, min_score=10):
    title = cfm_dict.pop("title")
    for spectra in cfm_dict:
        if not spectra is "fragments":
            data = cfm_dict[spectra]
            mass_intensities = data["mass_intensities"]
            spectra_data = spectra_dict[term_dict[spectra]]
            all_scores = data["scores"]
            for ms, peak_scores, spec, candidates in zip(mass_intensities, all_scores, spectra_data, data["candidates"]):
                if peak_scores:
                    best_index = np.argmax(peak_scores)
                    bestID = candidates[best_index]
                    if np.max(peak_scores) > min_score:
                        bar = plt.bar(ms[0], spec[1], color="r")[0]
                        bar_x = bar.get_x()
                        bar_w = bar.get_width()
                        bar_h = bar.get_height()
                        plt.text(bar_x + bar_w/2.,1.05*bar_h,str(bestID),ha='center', va='bottom')
                    else:
                        plt.bar(ms[0], spec[1])
                else:
                    plt.bar(ms[0], spec[1])

            plt.ylabel("intensity")
            plt.xlabel("mass/charge")
            plt.title(title+" - "+spectra)
            plt.savefig(output_path+"/cfm_spectra_"+spectra+".png")


def main():
    args = sys.argv
    if len(args) != 2:
        print "Usage: plot_spectra FILE_NAME"
        return

    data = io.read_spectra(args[1])
    plot_spectra(data)

if __name__ == '__main__':
    main()