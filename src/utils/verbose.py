import time
import sys

active = False
def initialise():
    global active
    global timestamps
    timestamps = {}
    active = True
    print "Verbose mode activated."

def info(name, *args):
    if active:
        if args:
            getattr(sys.modules[__name__],name)(*args)
        else:
            getattr(sys.modules[__name__],name)()

def show_setting(key, value):
    print "Set",key,"to:",value

def data_input(list_count,smiles_count):
        print "Read in",list_count,"list(s) with a total of",smiles_count,"SMILES strings."

def show_path_info(paths):
        distribution = [len(group) for group in paths]
        print "Found",sum(distribution),"unique paths."
        print "Path distribution:",distribution

def starting_length(length):
    print "Starting CS algorithm at length:",str(length)+"."

def current_length(length, name="CS"):
    print name,"search length:",str(length)+"."

def starting_struct_length(length):
    print "Starting Structure length:",str(length)+"."

def rep_struct_start():
    timestamps["rep_structs"] = time_in_ms()

def rep_struct_finish(structs):
    if len(structs):
        print "Found", len(structs),"rep. structures in",time_difference("rep_structs")+"."
        previous = -1
        buckets = []
        for struct in structs:
            rep_score = structs[struct]
            if previous == rep_score:
                buckets[-1][1] +=1
            else:
                buckets.append([rep_score, 0])
                previous = rep_score

        if buckets[-1][1] == 0:
            buckets.pop(-1)
        if buckets:
            print "Rep. tie-breaks: "+str(sum([bucket[1] for bucket in buckets]))

def rep_path_start():
    timestamps["rep_paths"] = time_in_ms()

def total_time_start():
    timestamps["start"] = time_in_ms()

def total_time_finish():
    print "Total calculation time:",str(time_difference("start"))+"."

def rep_path_finish(rep_path_count):
    if rep_path_count:
        print "Found", rep_path_count,"rep. paths in",time_difference("rep_paths")+"."

def start_paths():
    print "Started path search in molecules."
    timestamps["path_search"] = time_in_ms()

def finish_paths():
    print "Finished path search in",time_difference("path_search")+"."

def add_to_CS(name, times=0):
    if times <=1:
        print "Added path struct",name,"to CS."
    else:
        print "Added path struct",name,"to CS",times,"times."

def time_difference(term):
    return get_to_suitable_time(time_in_ms()-timestamps[term])

def time_in_ms():
    return int(round(time.time() * 1000))

#Assumes ms
def get_to_suitable_time(quantity):
    unit = "ms"

    if quantity > 1000:
        quantity /= 1000.0
        unit = "s"

        if quantity > 60:
            quantity /= 60.0
            unit = "min"

            if quantity > 60:
                quantity /= 60.0
                unit = "h"

    return str(quantity)+unit