from test_protected_sub6um_dc import run_both_states

if __name__ == '__main__':
    print("Testing t_pcm = 200 nm, L_c = 3.30 um, Gap = 98 nm...")
    run_both_states(resolution=15, l_c=3.30, gap=0.098, t_pcm=0.200)
