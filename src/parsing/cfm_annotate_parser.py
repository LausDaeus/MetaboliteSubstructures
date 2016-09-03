import src.utils.fileIO as io

def parse_cfm(data):
    parts = {}
    parts["title"] = data.pop(0).split()[-1]
    energy_type = ""
    index = 0
    for index,line in enumerate(data):
        if line in ['','\n','\r\n',' ']:
            break
        if any(char.isalpha() for char in line.split(" ")[0]) and len(line.split(" ")) == 1:
            energy_type = line.rstrip()
            print energy_type
            parts[energy_type] = {}
            parts[energy_type]["mass_intensities"] = []
            parts[energy_type]["candidates"] = []
            parts[energy_type]["scores"] = []
        else:
            peak_assignment = line.split(" ")
            parts[energy_type]["mass_intensities"].append([float(peak_assignment.pop(0)),float(peak_assignment.pop(0))])
            l = len(peak_assignment)
            candidates, scores = peak_assignment[:l/2], peak_assignment[l/2:]
            parts[energy_type]["candidates"].append([int(c) for c in candidates])
            parts[energy_type]["scores"].append([float(s.strip(")").strip("(")) for s in scores])

    index +=1

    if index < len(data):
        count = int(data[index])
        parts["fragments"] = []

        for i in xrange(index+1,count+index+1):
            parts["fragments"].append(data[i])

    return parts

