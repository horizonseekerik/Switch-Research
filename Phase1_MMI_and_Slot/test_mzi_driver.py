from sub3um_1x2_mzi_switch_3d import run_1x2_mzi_switch

if __name__ == '__main__':
    print("Testing 1x2 MZI Switch (3D FDTD, res=15)...")
    res_am = run_1x2_mzi_switch("amorphous", resolution=15)
    print("AMORPHOUS RESULT:")
    for k, v in res_am.items():
        print(f"  {k}: {v}")
    
    res_cr = run_1x2_mzi_switch("crystalline", resolution=15)
    print("CRYSTALLINE RESULT:")
    for k, v in res_cr.items():
        print(f"  {k}: {v}")
