import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn
import logging

from matplotlib import rcParams
from matplotlib.ticker import FormatStrFormatter


def ordering_for_types(D, types):
    possible_types = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    mutations = dict.fromkeys(possible_types)

    logging.debug("\n\nStarting mutations-order-based sorting of atoms..")
    for mut in mutations:
        mutations[mut] = np.where(types == mut)

    new_D = np.zeros_like(D)
    count = 0
    for mut in possible_types:
        logging.debug("Ordering type " + mut + " of mutation")
        new_D[:, count:count+len(mutations[mut][0])] = D[:, mutations[mut][0]]
        count += len(mutations[mut][0])
    return new_D


def remove_weak_mutations(mutational_catalogue, percentage):

    total_sum = np.sum(mutational_catalogue)

    cols_sum = np.sum(mutational_catalogue, axis=0)
    sorted_indices = np.argsort(cols_sum)

    partial_sum = 0
    for i in range(0, len(cols_sum)):
        partial_sum += cols_sum[sorted_indices[i]]
        if partial_sum/total_sum >= percentage:
            break

    cols_to_remove = np.sort(sorted_indices[0:i+1])
    new_mutational_catalogue = np.delete(mutational_catalogue,
                                         cols_to_remove, axis=1)

    return {"mutational_catalogue": new_mutational_catalogue,
            "removed_cols": cols_to_remove}


def add_removed_cols(dictionary, removed_cols):
    cols_dictionary = dictionary.shape[1]
    rows = dictionary.shape[0]
    cols = cols_dictionary + len(removed_cols)

    D = np.zeros((rows, cols))

    count = 0
    for i in range(0, cols):
        if count < len(removed_cols) and i == removed_cols[count]:
            count += 1
        else:
            index = i-count
            D[:, i] = np.real(dictionary[:, index])
    return D

def plot_genome_mutations2(dictionary):
    # get all mutations with context
    path = os.path.dirname(os.path.abspath(__file__))
    features = np.load("/home/veronica/Desktop/UVM/mutation_signatures/datasets/features.npy")
    # atoms normalization
    for r in range(dictionary.shape[0]):
        dictionary[r, :] /= np.sum(dictionary[r, :])

    cm = plt.cm.get_cmap('RdYlGn')
    colors = cm(np.linspace(0, 1, 6))

    _max = np.amax(np.amax(dictionary, axis=0)) + 0.05
    cols = 1
    rows = dictionary.shape[0]

    # start plot
    fig = plt.figure(figsize=(15,10))

    for k in range(dictionary.shape[0]):
        length = len(dictionary[k, :])
        x = np.asarray(range(0, length))
        w = dictionary[k, :]
        ax = fig.add_subplot(rows, cols, k + 1)

        if k == 0:  # saving first subplot for legend
            ax1 = ax

        weights = np.zeros(6)
        for j in range(6):
            i1 = 16 * j
            i2 = 16 * (j + 1)
            weights[j] = np.sum(w[i1:i2])

        N, bins, patches = plt.hist(x, bins=length, weights=w)
        width = patches[1].get_width()
        points = np.linspace(0, 95, 6)
        for k in range(6):
            ax.axvline(x=16*(k+1)*width, color='k', linestyle='dashed',
                       linewidth=1)

        N1, bins1, patches1 = plt.hist(points, bins=6, weights=weights,
                                       alpha=0.3)
        plt.xlim((0, dictionary.shape[1]))
        plt.ylim((0, _max))
        plt.xticks(bins+patches[0].get_width()/2 , features, rotation='vertical')
        rcParams.update({'font.size': 8})

        # colors labels
        labels = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
        for bin_size, _bin, patch in zip(N, bins, patches):
            for l in range(6):
                patches[16*(l+1)-1].set_label(labels[l])
                if _bin < 16 * (l+1)-1:
                    patch.set_facecolor(colors[l])
                    break

        for l in range(6):
            patches1[l].set_facecolor(colors[l])

    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='center',
               ncol=6,  mode='expand', borderaxespad=15.,
               numpoints=6, fontsize='xx-large')
    plt.tight_layout()
    plt.subplots_adjust(top=0.94)
    plt.show()


def plot_genome_mutations(dictionary):
    # get all mutations with context
    path = os.path.dirname(os.path.abspath(__file__))
    features = np.load("/home/veronica/Desktop/UVM/mutation_signatures/datasets/features.npy")
    # atoms normalization
    for r in range(dictionary.shape[0]):
        dictionary[r, :] /= np.sum(dictionary[r, :])

    cm = plt.cm.get_cmap('RdYlGn')
    colors = cm(np.linspace(0, 1, 6))

    _max = np.amax(np.amax(dictionary, axis=0)) + 0.05
    cols = 1
    rows = dictionary.shape[0]

    # start plot
    fig = plt.figure(figsize=(15,10))

    for k in range(dictionary.shape[0]):
        length = len(dictionary[k, :])
        x = np.asarray(range(0, length))
        w = dictionary[k, :]
        ax = fig.add_subplot(rows, cols, k + 1)

        if k == 0:  # saving first subplot for legend
            ax1 = ax

        weights = np.zeros(6)
        for j in range(6):
            i1 = 16 * j
            i2 = 16 * (j + 1)
            weights[j] = np.sum(w[i1:i2])

        N, bins, patches = plt.hist(x, bins=length, weights=w)
        width = patches[1].get_width()
        points = np.linspace(0, 95, 6)
        for k in range(6):
            ax.axvline(x=16*(k+1)*width, color='k', linestyle='dashed',
                       linewidth=1)

        N1, bins1, patches1 = plt.hist(points, bins=6, weights=weights,
                                       alpha=0.3)
        plt.xlim((0, dictionary.shape[1]))
        plt.ylim((0, _max))
        ax.xaxis.set_visible(False)
        ax.set_xticks(bins)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
        rcParams.update({'font.size': 5})

        # x labels
        threshold = (np.sum(np.array([rect.get_height()
                                      for rect in ax.patches]))
                     / len(ax.patches)) * 1.7
        for rect, label in zip(ax.patches, features):
            height = rect.get_height()
            if height > threshold:
                ax.text(rect.get_x() + rect.get_width() / 2, height + 0.02,
                        label, ha='center', va='bottom',
                        fontsize='large')

        # colors labels
        labels = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
        for bin_size, _bin, patch in zip(N, bins, patches):
            for l in range(6):
                patches[16*(l+1)-1].set_label(labels[l])
                if _bin < 16 * (l+1)-1:
                    patch.set_facecolor(colors[l])
                    break

        for l in range(6):
            patches1[l].set_facecolor(colors[l])

    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='center',
               ncol=6,  mode='expand', borderaxespad=15.,
               numpoints=6, fontsize='xx-large')
    plt.tight_layout()
    plt.subplots_adjust(top=0.94)
    plt.show()



def plot_coverage(dictionary, sample):
    # get all mutations with context
    features = np.load("/home/veronica/Desktop/UVM/mutation_signatures/datasets/features.npy")
    # atoms normalization
   # for r in range(dictionary.shape[0]):
     #   dictionary[r, :] /= np.sum(dictionary[r, :])

    cm = plt.cm.get_cmap('Spectral')
    colors = cm(np.linspace(0, 1, 6))

    _max = np.amax(sample) + 0.05
    cols = 1
    rows = dictionary.shape[0]

    # start plot
    fig = plt.figure(figsize=(15, 10))

    #for k in range(dictionary.shape[0]):
    length = len(dictionary)
    x = np.arange(0, length)
    w = dictionary
    #ax = fig.add_subplot(rows, cols, k + 1)

    #if k == 0:  # saving first subplot for legend
     #   ax1 = ax

    weights = np.zeros(6)
    for j in range(6):
        i1 = 16 * j
        i2 = 16 * (j + 1)
        weights[j] = np.sum(w[i1:i2])

    N, bins, patches = plt.hist(x, bins=length, weights=w)
    width = patches[1].get_width()
    points = np.linspace(0, 95, 6)
    for k in range(6):
        plt.axvline(x=16*(k+1)*width, color='k', linestyle='dashed',
                   linewidth=1)

    N1, bins1, patches1 = plt.hist(x, bins=length, weights=sample,
                                   alpha=0.5)
    plt.xlim((0, 96))
    plt.ylim((0, _max))
   # plt.xticks(bins+patches[0].get_width()/2 , features, rotation='vertical')
    rcParams.update({'font.size': 8})

    # colors labels
    labels = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    for bin_size, _bin, patch in zip(N, bins, patches):
        for l in range(6):
            patches[16*(l+1)-1].set_label(labels[l])
            if _bin < 16 * (l+1)-1:
                patch.set_facecolor(colors[l])
                break
    cm = plt.cm.get_cmap('bwr')
    colors = cm(np.linspace(0, 1, 6))
    for bin_size, _bin, patch in zip(N1, bins1, patches1):
        for l in range(6):
            if _bin < 16 * (l+1)-1:
                patch.set_facecolor(colors[l])
                break

    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='center',
               ncol=6,  mode='expand', borderaxespad=15.,
               numpoints=6, fontsize='xx-large')
    plt.tight_layout()
    plt.subplots_adjust(top=0.94)
    plt.show()


def plot_atom(dictionary):
    # get all mutations with context
    features = np.load("/home/veronica/Desktop/UVM/mutation_signatures/datasets/features.npy")
    # atoms normalization
   # for r in range(dictionary.shape[0]):
     #   dictionary[r, :] /= np.sum(dictionary[r, :])

    cm = plt.cm.get_cmap('RdYlGn')
    colors = cm(np.linspace(0, 1, 6))

    _max = 0.30#np.amax(dictionary) + 0.05
    cols = 1
    rows = dictionary.shape[0]

    # start plot
    fig = plt.figure(figsize=(15, 10))

    #for k in range(dictionary.shape[0]):
    length = len(dictionary)
    x = np.arange(0, length)
    w = dictionary
    #ax = fig.add_subplot(rows, cols, k + 1)

    #if k == 0:  # saving first subplot for legend
     #   ax1 = ax

    weights = np.zeros(6)
    for j in range(6):
        i1 = 16 * j
        i2 = 16 * (j + 1)
        weights[j] = np.sum(w[i1:i2])

    N, bins, patches = plt.hist(x, bins=length, weights=w)
    width = patches[1].get_width()
    points = np.linspace(0, 95, 6)
    for k in range(6):
        plt.axvline(x=16*(k+1)*width, color='k', linestyle='dashed',
                   linewidth=1)

    plt.xlim((0, 96))
    plt.ylim((0, _max))
    plt.xticks(bins+patches[0].get_width()/2 , features, rotation='vertical')
    rcParams.update({'font.size': 8})

    # colors labels
    labels = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    for bin_size, _bin, patch in zip(N, bins, patches):
        for l in range(6):
            patches[16*(l+1)-1].set_label(labels[l])
            if _bin < 16 * (l+1)-1:
                patch.set_facecolor(colors[l])
                break

    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='center',
               ncol=6,  mode='expand', borderaxespad=15.,
               numpoints=6, fontsize='xx-large')
    plt.tight_layout()
    plt.subplots_adjust(top=0.94)
    plt.show()
