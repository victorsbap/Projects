import math

print("MURO DE ARRIMO EM L")
print("\nVARIAVEIS DE ENTRADA")
gama_c = input("Entre gama c (kN/m³): ")
bm = input("Entre bm (m): ")
H = input("Entre H (m): ")
hb = input("Entre hb (m): ")
B = input("Entre B (m): ")
gama_s = input("Entre gama s (kN/m³): ")
phi = input("Entre angulo de atrito (graus): ")
SPT = input("Entre SPT medio: ")
Qsc_k = input("Entre Qsc k: ")

def peso_proprio(gama_c, bm, H, hb, B):
    print("\nPESO PROPRIO")
    gpar_k = float(gama_c) * float(bm)
    print("gpar_k =", gpar_k, "kN/m")
    Gm1_k = float(gpar_k) * float(H)
    print("Gm1_k =", Gm1_k, "kN")
    Gm2_k = float(gama_c) * float(bm) * float(hb)
    print("Gm2_k =", Gm2_k, "kN")
    gb_k = float(gama_c) * float(hb)
    print("gb_k =", gb_k, "kN/m")
    Gb_k = float(gama_c) * float(hb) * float (B)
    print("Gb_k =", Gb_k, "kN")
    return gpar_k, Gm1_k, gb_k, Gb_k, Gm2_k

def peso_solo(gama_s, B):
    print("\nPESO DO SOLO")
    sigma_vs = float(gama_s) * float(H)
    print("sigma_vs =", sigma_vs, "Pa")
    G_sk = float(sigma_vs) * float(B)
    print("G_sk =", G_sk, "kN")
    return sigma_vs, G_sk

def empuxo_solo(phi, gama_s, H, hb):
    print("\nEMPUXO DO SOLO")
    K = math.tan(math.radians(45 - (float(phi)/2)))**2
    print("K =", K)
    hs1 = float(K) * float(gama_s) * float(H)
    print("hs1 =", hs1, "kN/m²")
    hs2 = float(K) * float(gama_s) * (float(H)+float(hb))
    print("hs2 =", hs2, "kN/m²")
    Hs1_k = (float(hs1) * float(H))/2
    print("Hs1_k =", Hs1_k, "kN")
    Hs2_k = (float(K) * float(gama_s) * float(hb) * (2 * float(H) + float(hb))) / 2
    print("Hs2_k =", Hs2_k, "kN")
    Hs_k = float(Hs1_k) + float (Hs2_k)
    print("Hs_k =", Hs_k, "kN")
    return K, hs1, hs2, Hs1_k, Hs2_k, Hs_k

#Fatores de estabilidade
def deslizamento(phi, Gm1_k, Gm2_k, Gb_k, Gs_k, Hs_k):
    print("\nDESLIZAMENTO")
    mi_e = math.tan(math.radians(float(phi)*2/3))
    print("mi e =", mi_e)
    FS_d = (float(mi_e) * (float(Gm1_k)+float(Gm2_k)+float(Gb_k)+float(Gs_k))) / float(Hs_k)
    print("FS deslizamento =", FS_d)
    if FS_d >= 1.5:
        print("O MURO ESTA SEGURO QUANTO AO DESLIZAMENTO")
    else:
        print("O MURO NAO ESTA SEGURO QUANTO AO DESLIZAMENTO")
    return mi_e, FS_d

def tombamento(Gm1_k, Gm2_k, Gb_k, Gs_k, bm, B, Hs_k, H, hb):
    print("\nTOMBAMENTO")
    FS_t = ((float(Gm1_k) * float(Gm2_k)*(float(bm)/2)) + ((float(Gb_k) + float(Gs_k))*(float(bm)+(float(B)/2))))/(float(Hs_k)*((float(H)+float(hb))/3))
    print("FS tombamento =", FS_t)
    if FS_t >= 1.5:
        print("O MURO ESTA SEGURO QUANTO AO TOMBAMENTO")
    else:
        print("O MURO NAO ESTA SEGURO QUANTO AO TOMBAMENTO")
    return FS_t

def fundacoes(Hs_k, B, bm, Gm1_k, Gm2_k, Gb_k, Gs_k, H, hb, SPT):
    print("\nFUNDACOES")
    fat_k = float(Hs_k)/(float(B)+float(bm))
    print("fat k = ", fat_k, "kN/m")
    Nbase_k = float(Gm1_k) + float(Gm2_k) + float(Gb_k) + float(Gs_k)
    print("Nbase k = ", Nbase_k, "kN")
    Mbase_k = float(Hs_k) * (float(H)+float(hb))/3 + (float(Gm1_k)+float(Gm2_k))*(float(B)/2)- (float(Gb_k) + float(Gs_k)) * (float(hb)/2) #verificar se é hb ao inves de bm
    print("Mbase k = ", Mbase_k, "kN")
    A = float(B)+float(bm)
    print("Area da base = ", A, "m²")
    W = ((float(B)+float(bm))**2)/6
    print("W = ", W, "m³")
    e = float(Mbase_k)/float(Nbase_k)
    print("e = ", e, "m")
    #e_lim = (float(B)+float(bm))/6
    e_lim = float(A) / 6
    print("e limite = ", e_lim, "m")
    sigma_sk_max = (float(Nbase_k)/float(A)) + (float(Mbase_k)/float(W))
    print("sigma sk maximo = ", sigma_sk_max, "kN/m²")
    sigma_sk_min = (float(Nbase_k)/float(A)) - (float(Mbase_k)/float(W))
    print("sigma sk minimo = ", sigma_sk_min, "kN/m²")
    sigma_s_adm = float(SPT)/50
    print("sigma s admissivel = ", sigma_s_adm, "MPa")
    if e <= e_lim:
        print("\nNAO HAVERA LEVANTAMENTO DE BASE")
    else:
        print("\nHAVERA LEVANTAMENTO DE BASE")
    if (float(sigma_sk_max)/1000) <= float(sigma_s_adm):
        print("E O SOLO TEM CAPACIDADE DE SUPORTE")
    else:
        print("E O SOLO NAO POSSUI CAPACIDADE DE SUPORTE")
    return fat_k, Nbase_k, Mbase_k, A, W, sigma_sk_max, sigma_sk_min, sigma_s_adm, e, e_lim

def esforcos(Gm1_k, Hs2_k, B, Hs_k, bm, sigma_sk_min, sigma_sk_max, Hs1_k, Gb_k, hb, Gs_k, Qsc_k, hs1, hs2):
    fat_k = float(Hs_k) / (float(B) + float(bm))
    print("\nFORÇA NORMAL")
    Na_k = 0
    print("NA k = ", Na_k, "kN")
    Nb_k = -float(Gm1_k)
    print("NB k = ", Nb_k, "kN")
    Nc_k = (-float(Hs2_k) + (float(B)*float(fat_k)))
    print("NC k = ", Nc_k, "kN")
    Nd_k = (-float(Hs2_k))
    print("ND k = ", Nd_k, "kN")
    print("\nFORÇA CORTANTE")
    Va_k = 0
    print("VA = ", Va_k, "kN")
    Vb_k = -float(Hs1_k)
    print("VB k = ", Vb_k, "kN")
    sigma_sp = float(sigma_sk_min) + (float(B)*((float(sigma_sk_max)-float(sigma_sk_min))))/(float(B)+float(bm))
    print("sigma sp = ", sigma_sp, "kN/m²")
    R_sp = float(B)*((float(sigma_sp)+float(sigma_sk_min)))/2
    print("R sp = ", R_sp, "kN")
    Vc_k = float(Gb_k) + float(Gs_k) - float(R_sp)
    print("VC k = ", Vc_k, "kN")
    Vd_k = 0
    print("VD k = ", Vd_k)
    print("\nMOMENTO") #Última vez alterado
    Ma_k = 0
    print("MA k = ", Ma_k, "kN.m")
    Mb_k = (float(Hs1_k)*(float(H)))/3
    print("MB k = ", Mb_k, "kN.m")
    Mc_k = (float(R_sp) * float(B) * (float(sigma_sp) + (2 * float(sigma_sk_min))) / (3 * (float(sigma_sp) + float(sigma_sk_min)))) + (((float(fat_k) * float(B) * float(hb))) / 2) - ((float(Gb_k) + float(Gs_k) + float(Qsc_k)) * ((float(B) / 2))) - (float(Hs2_k) * ((((float(hb) / 2))) - (float(hb) * (float(hs2) + (float(hs1) * 2)) / (3 * (float(hs2) + float(hs1))))))
    print("MC k = ", Mc_k, "kN.m")
    Md_k = -float(Hs2_k)*((float(hb)/2)-(float(hb)*(float(hs2)+(2*float(hs1))))/(3*(float(hs2)+float(hs1))))
    print("MD k = ", Md_k, "kN.m")
    return Na_k, Nb_k, Nc_k, Nd_k, Va_k, Vb_k, sigma_sp, R_sp, Vc_k, Vd_k, Ma_k, Mb_k, Mc_k, Md_k

gpar_k, Gm1_k, gb_k, Gb_k, Gm2_k = peso_proprio(gama_c, bm, H, hb, B)
sigma_vs, Gs_k = peso_solo(gama_s, B)
K, hs1, hs2, Hs1_k, Hs2_k, Hs_k = empuxo_solo(phi, gama_s, H, hb)
fat_k, Nbase_k, Mbase_k, A, W, sigma_sk_max, sigma_sk_min, sigma_s_adm, e, e_lim = fundacoes(Hs_k, B, bm, Gm1_k, Gm2_k, Gb_k, Gs_k, H, hb, SPT)
#Na_k, Nb_k, Nc_k, Nd_k, Va_k, Vb_k, Vc_k, Vd_k, sigma_sp, R_sp = esforcos(Gm1_k, Hs2_k, B, Hs_k, bm, sigma_sk_min, sigma_sk_max) #Última modificação
estabilidade_desl = deslizamento(phi, Gm1_k, Gm2_k, Gb_k, Gs_k, Hs_k)
estabilidade_tomb = tombamento(Gm1_k, Gm2_k, Gb_k, Gs_k, bm, B, Hs_k, H, hb)
acoes_fundacao = fundacoes(Hs_k, B, bm, Gm1_k, Gm2_k, Gb_k, Gs_k, H, hb, SPT)
esforcos_internos = esforcos(Gm1_k, Hs2_k, B, Hs_k, bm, sigma_sk_min, sigma_sk_max, Hs1_k, Gb_k, hb, Gs_k, Qsc_k, hs1, hs2)

result_pp = peso_proprio(gama_c, bm, H, hb, B)
result_ps = peso_solo(gama_s, B)
result_empuxo = empuxo_solo(phi, gama_s, H, hb)
estabilidade_desl = deslizamento(phi, Gm1_k, Gm2_k, Gb_k, Gs_k, Hs_k)
estabilidade_tomb = tombamento(Gm1_k, Gm2_k, Gb_k, Gs_k, bm, B, Hs_k, H, hb)
acoes_fundacao = fundacoes(Hs_k, B, bm, Gm1_k, Gm2_k, Gb_k, Gs_k, H, hb, SPT)
esforcos_internos = esforcos(Gm1_k, Hs2_k, B, Hs_k, bm, sigma_sk_min, sigma_sk_max, Hs1_k, Gb_k, hb, Gs_k, Qsc_k, hs1, hs2)