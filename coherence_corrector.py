
#If the coherence is above 1 or below 0, it is set to 1 or 0 respectively.
def coherence_corrector():
    coherence=np.array(modulus_G_sqrd/(np.array(cs_real_mean)*np.array(cs_ith_real_mean_array)))
        print('coherence',coherence)

        new_coherence = []
        for i in coherence:
            if i > 1:
                new_coherence.append(1)
            elif i < 0:
                new_coherence.append(0)
            else:
                new_coherence.append(i)
    
        new_coherence=np.array(new_coherence)