import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

if __name__ == '__main__':

    def genetic_drift(n=1000, freq=0.5, g=1000):
        data = [[], []]
        p_A = freq
        fixation_gen = 0
        for i in range(g):
            fixation_gen = i
            data[ 0 ].append(i)
            data[ 1 ].append(p_A)
            p_A = np.random.binomial(n,
                                     p_A) / n  #
            if p_A == 0 or p_A == 1:
                break
        return data, fixation_gen


    def genetic_drift_with_selection(n=100000, freq=0.5, h=0.25, s=0.2, g=100):
        data = [ ]
        genotype_data = [ ]
        p_A = freq
        p_a = 1 - freq
        w_AA = 1
        w_Aa = 1 - (h * s)
        w_aa = 1 - s
        #Hardy Weinberg equilibrium
        p_AA = p_A * p_A
        p_Aa = 2 * p_A * p_a
        p_aa = p_a * p_a
        fixation_gen = 0
        for i in range(g):  # Repeat for all generations
            wbar = (p_AA * w_AA) + (p_Aa * w_Aa) + (p_aa * w_aa)
            # New frequencies after selection
            p_AA_wAA = p_AA * w_AA / wbar
            p_Aa_wAa = p_Aa * w_Aa / wbar
            p_aa_waa = p_aa * w_aa / wbar  # Eq. to  1 - p_AA_wAA - p_Aa_wAa
            # Draw one sample of the population following a multinomial distribution
            new_population = np.random.multinomial(n, [ p_AA_wAA, p_Aa_wAa, p_aa_waa ], size=1)
            genotype_data.append(new_population[ 0 ] / n)
            p_AA = new_population[ 0 ][ 0 ] / n
            p_Aa = new_population[ 0 ][ 1 ] / n
            p_aa = new_population[ 0 ][ 2 ] / n
            p_A = p_AA + 0.5 * p_Aa  # Calculate the new allele frequencies
            p_a = p_aa + 0.5 * p_Aa
            data.append([ p_A ])
            fixation_gen = i
            if p_A == 0 or p_A == 1:  # Finish when the allele reaches fixation
                break
        return data, genotype_data, fixation_gen

    p = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('-n',
        dest = 'population_size',
        help = 'Population size',
        type = int,
        default = 100)

    p.add_argument('-f',
        dest = 'freq',
        help = 'Frequency of allele',
        type = float,
        default = 0.5)

    p.add_argument('-g',
        dest = 'generations',
        help = 'Maximum number of generations to simulate',
        type = int,
        default = 1000)

    p.add_argument('-r',
        dest = 'repetitions',
        help = 'Number of simulations',
        type = int,
        default = 3)

    p.add_argument('--with_selection', action="store_true", help = 'Add selection to the model')
    p.add_argument('-s',
    dest = 'selection_coef',
    help = 'Selection coefficient (When the model applies selection)',
    type = float,
    default = 0.2)

    p.add_argument('-d',
    dest = 'dominance_coef',
    help = 'Dominance coefficient (h) (When the model applies selection)',
    type = float,
    default = 0.25)

    args = p.parse_args()

    if not args.with_selection: # Simulating genetic drift without selection
        # Get parameters
        rep = args.repetitions
        n = args.population_size
        freq = args.freq
        g = args.generations
        pop_freq = []
        fixation = []

        # Initialize plot
        plt.figure(figsize=(12, 10))
        ax1 = plt.subplot(1,1,1)
        cmap = plt.get_cmap('viridis')
        for i in range(rep): # Simulate for each repetition
            sim, fixation_gen = genetic_drift(n=n, freq=0.5, g=g)
            pop_freq.append(sim[1][-1]) #frequency and fixation generation to an array to show the average
            if fixation_gen != g - 1:
                fixation.append(fixation_gen)

            ax1.plot(sim[0], sim[1], c=cmap(float(i)/rep))
        #Setting labels and titles
        ax1.set_title('Genetic drift: Frequency of allele over up to {} generations\n (N={}, freq in gen 0={}, {} '
                      'repetitions)'.format(g, n, freq, rep))
        ax1.set_ylabel('Frequency')
        ax1.set_xlabel('Generation')
        ax1.set_ylim([0, 1])
        f_text = 'Average frequency in last generation: {}'.format(round(sum(pop_freq)/len(pop_freq),3))
        g_text = 'Average number of generations until fixation of a gene: {}'.format(int(sum(fixation)/len(fixation)) if len(fixation) > 0 else 'None of the repetitions reached fixation')
        ax1.text(0.5,-0.1, f'{f_text}\n{g_text}', size=12, ha="center", transform=ax1.transAxes)
        plt.show()

    if args.with_selection:
        rep = args.repetitions
        n = args.population_size
        freq = args.freq
        h = args.dominance_coef
        s = args.selection_coef
        g = args.generations
        pop_freq = []
        fixation = []

        fig = plt.figure(figsize=(12, 10))
        ax1 = plt.subplot(2,1,1)
        ax1.set_title('Freq of allele in a Wright-Fisher model with selection\n (N={}, h={}, s = {}, freq in gen 0={}, {} repetitions)'.format(n, h, s, freq, rep))
        ax1.set_ylabel('Frequency')
        ax1.set_xlabel('Generation')
        cmap = plt.get_cmap('viridis') # Apply the color palette
        for i in range(rep): #Simulate for each repetition
            data, genotype_data, fixation_gen = genetic_drift_with_selection(n = n, freq = freq, h = h, s = s, g = g)
            pop_freq.append(data[-1][-1]) #last frequency and fixation generation to an array to show the average
            if fixation_gen != g - 1:
                fixation.append(fixation_gen)
            ax1.plot(data, c=cmap(float(i)/rep))

        ax2 = plt.subplot(2,1,2) #subplot for the genotype frequency
        p1 = plt.bar([], [])
        rects = plt.bar(['AA', 'Aa', 'aa'], genotype_data[0], color=cmap(float(i)/rep))
        ax2.set_ylim(0, 1)
        ax2.set_title('Frequency of genotype')
        ax2.set_ylabel('Frequency')
        ax2.set_xlabel('Genotype')


        def animate(i): # Animate function, to update the bar plot
            for barchart, yi in zip(rects, genotype_data[i]):
                barchart.set_height(yi)
                ax2.set_title('Frequency of genotypes in generation {} (Last repetition)'.format(i))

            return rects

        f_text = 'Average frequency in last generation: {}'.format(round(sum(pop_freq)/len(pop_freq),3))
        g_text = 'Average number of generations until fixation of a mutant gene'.format(int(sum(fixation)/len(fixation))
                                                                                        if len(fixation) > 0 else 'None of the repetitions reached fixation')
        ax2.text(0.5,-0.23, f'{f_text}\n{g_text}', size=12, ha="center", transform=ax2.transAxes)
        anim = animation.FuncAnimation(fig, animate, frames=len(genotype_data), interval=10, repeat=False)
        plt.show()


