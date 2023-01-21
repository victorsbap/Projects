import math

#variáveis de entrada
gama_c = input("Entre gama c: ")
bm = input("Entre bm: ")
H = input("Entre H: ")
hb = input("Entre hb: ")
B = input("Entre B: ")
gama_s = input("Entre gama s: ")
phi = input("Entre angulo de atrito: ")
SPT = input("Entre SPT medio: ")

def peso_proprio(gama_c, bm, H, hb, B):
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
    sigma_vs = float(gama_s) * float(H)
    print("sigma_vs =", sigma_vs, "Pa")
    G_sk = float(sigma_vs) * float(B)
    print("G_sk =", G_sk, "kN")
    return sigma_vs, G_sk

def empuxo_solo(phi, gama_s, H, hb):
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
    FS_t = ((float(Gm1_k) * float(Gm2_k)*(float(bm)/2)) + ((float(Gb_k) + float(Gs_k))*(float(bm)+(float(B)/2))))/(float(Hs_k)*((float(H)+float(hb))/3))
    print("FS tombamento =", FS_t)
    if FS_t >= 1.5:
        print("O MURO ESTA SEGURO QUANTO AO TOMBAMENTO")
    else:
        print("O MURO NAO ESTA SEGURO QUANTO AO TOMBAMENTO")
    return FS_t

def fundacoes(Hs_k, B, bm, Gm1_k, Gm2_k, Gb_k, Gs_k, H, hb, SPT):
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
        print("NAO HAVERA LEVANTAMENTO DE BASE")
    else:
        print("HAVERA LEVANTAMENTO DE BASE")
    if (float(sigma_sk_max)/1000) <= float(sigma_s_adm):
        print("O SOLO TEM CAPACIDADE DE SUPORTE")
    else:
        print("E O SOLO NAO POSSUI CAPACIDADE DE SUPORTE")
    return fat_k, Nbase_k, Mbase_k, A, W, sigma_sk_max, sigma_sk_min, sigma_s_adm, e, e_lim

#def esforços():
    #Na_k = 0
    #print("Na k = ", Na_k)
    #Nb_k = -float(Gm1_k)

gpar_k, Gm1_k, gb_k, Gb_k, Gm2_k = peso_proprio(gama_c, bm, H, hb, B)
sigma_vs, Gs_k = peso_solo(gama_s, B)
K, hs1, hs2, Hs1_k, Hs2_k, Hs_k = empuxo_solo(phi, gama_s, H, hb)
estabilidade_desl = deslizamento(phi, Gm1_k, Gm2_k, Gb_k, Gs_k, Hs_k)
estabilidade_tomb = tombamento(Gm1_k, Gm2_k, Gb_k, Gs_k, bm, B, Hs_k, H, hb)
acoes_fundacao = fundacoes(Hs_k, B, bm, Gm1_k, Gm2_k, Gb_k, Gs_k, H, hb, SPT)

result_pp = peso_proprio(gama_c, bm, H, hb, B)
result_ps = peso_solo(gama_s, B)
result_empuxo = empuxo_solo(phi, gama_s, H, hb)
estabilidade_desl = deslizamento(phi, Gm1_k, Gm2_k, Gb_k, Gs_k, Hs_k)
estabilidade_tomb = tombamento(Gm1_k, Gm2_k, Gb_k, Gs_k, bm, B, Hs_k, H, hb)
acoes_fundacao = fundacoes(Hs_k, B, bm, Gm1_k, Gm2_k, Gb_k, Gs_k, H, hb, SPT)