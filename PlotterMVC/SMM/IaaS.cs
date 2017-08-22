using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace PlotterMVC.SMM
{
    public class IaaS
    {
        public double t_min { get; private set; }
        public double t_max { get; set; }
        double t_step { get; set; }

        public double gamma_min { get; set; }
        public double gamma_max { get; private set; }
        double gamma_step { get; set; }

        public double mu { get; private set; }

        double mu1;
        double mu2;
        double tcr;
        double k;
        double p75;

        int n;
        int m;

        int n_max;
        int m_max;

        public double[] gamma { get; private set; }
        public double[] T { get; private set; }

        public double[,] P0 { get; private set; }
        public double[,] P1 { get; private set; }
        public double[,] P2 { get; private set; }
        public double[,] P3 { get; private set; }
        public double[,] P4 { get; private set; }
        public double[,] P5 { get; private set; }
        public double[,] P6 { get; private set; }
        public double[,] P7 { get; private set; }
        public double[,] P8 { get; private set; }
        public double[,] P9 { get; private set; }
        public double[,] P10 { get; private set; }
        public double[,] P11 { get; private set; }
        public double[,] P12 { get; private set; }
        public double[,] P13 { get; private set; }

        public double[,] AVL { get; private set; }
        public double[,] UNAVL { get; private set; }

        public IaaS()
        {
            t_min = 1;
            t_max = 280;
            t_step = 1;

            gamma_min = 0.000005;
            gamma_max = 0.00001;
            gamma_step = 0.000001;

            n_max = 280;
            m_max = 100;

            mu = 8;
            tcr = 0.1666;
            k = 1;
            p75 = 1;            
        }

        private void Init()
        {
            n = Calc_n(t_min, t_max, t_step);
            m = Calc_m(gamma_min, gamma_max, gamma_step);

            gamma = new double[m];
            T = new double[n];

            P0 = new double[m, n];
            P1 = new double[m, n];
            P2 = new double[m, n];
            P3 = new double[m, n];
            P4 = new double[m, n];
            P5 = new double[m, n];
            P6 = new double[m, n];
            P7 = new double[m, n];
            P8 = new double[m, n];
            P9 = new double[m, n];
            P10 = new double[m, n];
            P11 = new double[m, n];
            P12 = new double[m, n];
            P13 = new double[m, n];

            AVL = new double[m, n];
            UNAVL = new double[m, n];
        }
        
        public void Set_T_min(double Tmin)
        {
            this.t_min = Tmin;
        }
        
        public void Set_T_max(double Tmax)
        {
            this.t_max = Tmax;
        }

        public void Set_gamma_min(double gamma_min)
        {
            this.gamma_min = gamma_min;
        }

        public void Set_gamma_max(double gamma_max)
        {
            this.gamma_max = gamma_max;
        }
        
        public void Set_mu(double mu)
        {
            this.mu = mu;
        }

        public void Calculate()
        {
            Init();

            mu1 = Calc_mu1(mu);
            mu2 = Calc_mu2(mu);

            for (int i = 0; i < m; i++)
            {
                gamma[i] = Calc_gamma(gamma_max, gamma_step, i);
                double gamma1 = Calc_gamma1(k, gamma[i]); // Failure rate of hot and warm PMs
                double gamma2 = Calc_gamma2(gamma1); // Hidden failure rate of hot and warm PMs before

                // second control of technical state (TSC)
                double gamma3 = Calc_gamma3(gamma2); // Hidden failure rate of hot and warm PMs after TSC

                for (int j = 0; j < n; j++)
                {
                    T[j] = Calc_T(t_min, t_step, j); // Time of use IaaS Cloud

                    //for ought state
                    double t0 = Calc_t0(gamma2, T[j]); // Duration of time

                    //for first state
                    double t1 = Calc_t1(gamma1, tcr); // Duration of time

                    //for second state
                    double t2 = Calc_t2(gamma2, mu); // Duration of time

                    //for third state
                    double t3 = Calc_t3(t2); // Duration of time

                    //for fourth state
                    double t4 = Calc_t4(gamma1, tcr); // Duration of time

                    //for fifth state
                    double t5 = Calc_t5(mu2); // Duration of time

                    //for sixth state
                    double t6 = Calc_t6(t4); // Duration of time

                    //for seventh state
                    double t7 = Calc_t7(mu1); // Duration of time

                    //Inter mediate equation for determening of eighth time duration
                    double t08 = Calc_t08(T[j], gamma2);

                    //for eighth state
                    double t8 = Calc_t8(T[j], t08); // Duration of time

                    //for ninth state
                    double t9 = Calc_t9(gamma3, tcr); // Duration of time

                    //for tenth state
                    double t10 = Calc_t10(t8); // Duration of time

                    //for eleventh state
                    double t11 = Calc_t11(gamma3, tcr); // Duration of time

                    //for twelfth state
                    double t12 = Calc_t12(t8); // Duration of time

                    //for thirteenth stat
                    double t13 = Calc_t13(t11); // Duration of time   

                    // ======================== Transient probabilities ========================
                    double p01 = Calc_p01(gamma2, T[j]);
                    double p08 = Calc_p08(p01);
                    double p10 = Calc_p10(gamma1, tcr);
                    double p12 = Calc_p12(gamma1, T[j]);
                    double p13 = Calc_p13(p10, p12);
                    double p92 = Calc_p92(gamma3, T[j]);
                    double p98 = Calc_p98(gamma3, tcr);
                    double p93 = Calc_p93(p92, p98);
                    double p24 = Calc_p24(gamma2, mu1, T[j]);
                    double p34 = Calc_p34(p24);
                    double p210 = Calc_p210(gamma2, mu1, T[j]);
                    double p312 = Calc_p312(p210);
                    double p20 = Calc_p20(p24, p210);
                    double p30 = Calc_p30(p20);
                    double p42 = Calc_p42(gamma1, tcr);
                    double p43 = Calc_p43(p42);
                    double p65 = Calc_p65(p42);
                    double p45 = Calc_p45(p42, p43);
                    double E1 = Calc_E1(mu2, T[j]);
                    double E2 = Calc_E2(mu2, T[j], E1);
                    double p52 = Calc_p52(mu2, T[j], E1, E2);
                    double p56 = Calc_p56(mu2, T[j]);
                    double p53 = Calc_p53(p52, p56);
                    double p67 = Calc_p67(p65);
                    double p1110 = Calc_p1110(gamma3, tcr);
                    double p115 = Calc_p115(p1110);
                    double p1312 = Calc_p1312(p1110);
                    double p135 = Calc_p135(p115);

                    // ======================== Intermediate parameters ========================
                    double a1 = Calc_a1(p10);
                    double a2 = Calc_a2(p20);
                    double a3 = Calc_a3(p30);
                    double a4 = Calc_a4(p01);
                    double a5 = Calc_a5(p12);
                    double a6 = Calc_a6(p42);
                    double a7 = Calc_a7(p52);
                    double a8 = Calc_a8(p92);
                    double a9 = Calc_a9(p13);
                    double a10 = Calc_a10(p43);
                    double a11 = Calc_a11(p53);
                    double a12 = Calc_a12(p93);
                    double a13 = Calc_a13(p24);
                    double a14 = Calc_a14(p34);
                    double a15 = Calc_a15(p45);
                    double a16 = Calc_a16(p65);
                    double a17 = Calc_a17(p75);
                    double a18 = Calc_a18(p115);
                    double a19 = Calc_a19(p135);
                    double a20 = Calc_a20(p56);
                    double a21 = Calc_a21(p67);
                    double a22 = Calc_a22(p08);
                    double a23 = Calc_a23(p98);
                    double a24 = Calc_a24(p210);
                    double a25 = Calc_a25(p1110);
                    double a26 = Calc_a26(p312);
                    double a27 = Calc_a27(p1312);
                    double k7 = Calc_k7(a20, a21);
                    double k8 = Calc_k8(a22, a23);
                    double k10 = Calc_k10(a24, a25);
                    double k12 = Calc_k12(a26, a27);
                    double k20 = Calc_k20(a4, a5);
                    double k21 = Calc_k21(a6, a13);
                    double k22 = Calc_k22(a6, a14);
                    double k23 = Calc_k23(a7, a15);
                    double k28 = Calc_k28(a8, k8);
                    double k24 = Calc_k24(a7, a16, a20);
                    double k25 = Calc_k25(a7, a17, k7);
                    double k26 = Calc_k26(a7, a18, k10);
                    double k27 = Calc_k27(a7, a19, k12);
                    double k29 = Calc_k29(a13, k23);
                    double k30 = Calc_k30(a14, k23);
                    double k31 = Calc_k31(a4, a9);
                    double k32 = Calc_k32(a10, a13);
                    double k33 = Calc_k33(a10, a14);
                    double k34 = Calc_k34(a11, a15);
                    double k35 = Calc_k35(a11, a16, a20);
                    double k36 = Calc_k36(a11, a17, k7);
                    double k37 = Calc_k37(a11, a18, k10);
                    double k38 = Calc_k38(a11, a19, k12);
                    double k39 = Calc_k39(a12, k8);
                    double k50 = Calc_k50(a16, a20);
                    double k51 = Calc_k51(a17, k7);
                    double k52 = Calc_k52(a18, k10);
                    double k53 = Calc_k53(a19, k12);
                    double k54 = Calc_k54(a13, a15);
                    double k55 = Calc_k55(a14, a15);
                    double Zn5 = Calc_Zn5(k50, k51);
                    double c20 = Calc_c20(k24, k52);
                    double c21 = Calc_c21(k24, k54);
                    double c22 = Calc_c22(k24, k53);
                    double c23 = Calc_c23(k24, k55);
                    double c24 = Calc_c24(k25, k52);
                    double c25 = Calc_c25(k25, k54);
                    double c26 = Calc_c26(k25, k53);
                    double c27 = Calc_c27(k25, k55);
                    double c28 = Calc_c28(k20, Zn5);
                    double c29 = Calc_c29(k28, Zn5);
                    double c30 = Calc_c30(k22, Zn5);
                    double c31 = Calc_c31(k27, Zn5);
                    double c32 = Calc_c32(k30, Zn5);
                    double c33 = Calc_c33(a13, c28);
                    double c34 = Calc_c34(a13, c29);
                    double c35 = Calc_c35(a13, c22);
                    double c36 = Calc_c36(a13, c23);
                    double c37 = Calc_c37(a13, c26);
                    double c38 = Calc_c38(a13, c27);
                    double c39 = Calc_c39(a13, c30);
                    double c40 = Calc_c40(a13, c31);
                    double c41 = Calc_c41(a13, c32);
                    double Zn2 = Calc_Zn2(Zn5, c20, c21, c24, c25, k21, k26, k29);
                    double c42 = Calc_c42(a14, Zn2);
                    double c43 = Calc_c43(k52, c28);
                    double c44 = Calc_c44(k52, c29);
                    double c45 = Calc_c45(k52, c22);
                    double c46 = Calc_c46(k52, c23);
                    double c47 = Calc_c47(k52, c26);
                    double c48 = Calc_c48(k52, c27);
                    double c49 = Calc_c49(k52, c30);
                    double c50 = Calc_c50(k52, c31);
                    double c51 = Calc_c51(k52, c32);
                    double c52 = Calc_c52(k54, c28);
                    double c53 = Calc_c53(k54, c29);
                    double c54 = Calc_c54(k54, c22);
                    double c55 = Calc_c55(k54, c23);
                    double c56 = Calc_c56(k54, c26);
                    double c57 = Calc_c57(k54, c27);
                    double c58 = Calc_c58(k54, c30);
                    double c59 = Calc_c59(k54, c31);
                    double c60 = Calc_c60(k54, c32);
                    double c61 = Calc_c61(k53, Zn2);
                    double c62 = Calc_c62(k55, Zn2);
                    double kf41 = Calc_kf41(c33, c34);
                    double kf42 = Calc_kf42(c35, c36, c37, c38, c39, c40, c41, c42);
                    double kf51 = Calc_kf51(c43, c44, c52, c53);
                    double kf52 = Calc_kf52(c45, c46, c47, c48, c49, c50, c51, c54, c55, c56, c57, c58, c59, c60, c61, c62);
                    double kf300 = Calc_kf300(Zn2, Zn5);
                    double kf301 = Calc_kf301(k31, Zn2, Zn5);
                    double kf302 = Calc_kf302(k39, Zn2, Zn5);
                    double kf303 = Calc_kf303(k32, c28, Zn5);
                    double kf304 = Calc_kf304(k32, c29, Zn5);
                    double kf305 = Calc_kf305(k32, c22, Zn5);
                    double kf306 = Calc_kf306(k32, c23, Zn5);
                    double kf307 = Calc_kf307(k32, c26, Zn5);
                    double kf308 = Calc_kf308(k32, c27, Zn5);
                    double kf309 = Calc_kf309(k32, c30, Zn5);
                    double kf310 = Calc_kf310(k32, c31, Zn5);
                    double kf311 = Calc_kf311(k32, c32, Zn5);
                    double kf312 = Calc_kf312(k37, c28, Zn5);
                    double kf313 = Calc_kf313(k37, c29, Zn5);
                    double kf314 = Calc_kf314(k37, c22, Zn5);
                    double kf315 = Calc_kf315(k37, c23, Zn5);
                    double kf316 = Calc_kf316(k37, c26, Zn5);
                    double kf317 = Calc_kf317(k37, c27, Zn5);
                    double kf318 = Calc_kf318(k37, c30, Zn5);
                    double kf319 = Calc_kf319(k37, c31, Zn5);
                    double kf320 = Calc_kf320(k37, c32, Zn5);
                    double kf321 = Calc_kf321(k33, Zn2, Zn5);
                    double kf322 = Calc_kf322(k38, Zn2, Zn5);
                    double kf323 = Calc_kf323(k34, kf41, Zn5);
                    double kf324 = Calc_kf324(k34, kf42, Zn5);
                    double kf325 = Calc_kf325(k35, kf51);
                    double kf326 = Calc_kf326(k35, kf52);
                    double kf327 = Calc_kf327(kf325);
                    double kf328 = Calc_kf328(k36, kf52);
                    double Ch3 = Calc_Ch3(kf301, kf302, kf303, kf304, kf312, kf313, kf323, kf325, kf327);
                    double Zn3 = Calc_Zn3(kf300, kf305, kf306, kf307, kf308, kf309, kf310, kf311, kf314, kf315, kf316, kf317, kf318, kf319, kf320, kf321, kf322, kf324, kf326, kf328);
                    double kf200 = Calc_kf200(Zn2, Zn3);
                    double kf201 = Calc_kf201(c28, Zn3);
                    double kf202 = Calc_kf202(c29, Zn3);
                    double kf203 = Calc_kf203(c22, Ch3);
                    double kf204 = Calc_kf204(c23, Ch3);
                    double kf205 = Calc_kf205(c26, Ch3);
                    double kf206 = Calc_kf206(c27, Ch3);
                    double kf207 = Calc_kf207(c30, Ch3);
                    double kf208 = Calc_kf208(c31, Ch3);
                    double kf209 = Calc_kf209(c32, Ch3);
                    double Ch2 = Calc_Ch2(kf201, kf202, kf203, kf204, kf205, kf206, kf207, kf208, kf209);
                    double kf100 = Calc_kf100(k10, Ch2, kf200);
                    double kf120 = Calc_kf120(k12, kf300);
                    double kf400 = Calc_kf400(a13, Ch2, Zn3);
                    double kf401 = Calc_kf401(a14, Ch3, kf200);
                    double kf402 = Calc_kf402(kf400, kf401, kf200, Zn3);
                    double kf500 = Calc_kf500(kf51, Zn3);
                    double kf501 = Calc_kf501(kf52, Ch3);
                    double kf502 = Calc_kf502(Zn2, Zn3, Zn5);
                    double kf503 = Calc_kf503(kf500, kf501, kf502);
                    double kf600 = Calc_kf600(a20, kf503);
                    double kf700 = Calc_kf700(a21, kf600);
                    double kf800 = Calc_kf800(a22, k8);
                    double kf201_1 = Calc_kf201_1(Ch2, kf200);
                    double U = Calc_U(a4, kf201_1, kf300, kf402, kf503, kf600, kf700, kf800, kf100, kf120);

                    //======================= Additional steady state probabilities =============
                    double Pi0 = Calc_Pi0(U);
                    double Pi1 = Calc_Pi1(a4, Pi0);
                    double Pi2 = Calc_Pi2(kf201_1, Pi0);
                    double Pi3 = Calc_Pi3(kf300, Pi0);
                    double Pi4 = Calc_Pi4(kf402, Pi0);
                    double Pi5 = Calc_Pi5(kf503, Pi0);
                    double Pi6 = Calc_Pi6(kf600, Pi0);
                    double Pi7 = Calc_Pi7(kf700, Pi0);
                    double Pi8 = Calc_Pi8(kf800, Pi0);
                    double Pi9 = Calc_Pi9(Pi8);
                    double Pi10 = Calc_Pi10(kf100, Pi0);
                    double Pi11 = Calc_Pi11(Pi10);
                    double Pi12 = Calc_Pi12(kf120, Pi0);
                    double Pi13 = Calc_Pi13(Pi12);

                    double A0 = Calc_A0(Pi0, t0);
                    double A1 = Calc_A1(Pi1, t1);
                    double A2 = Calc_A2(Pi2, t2);
                    double A3 = Calc_A3(Pi3, t3);
                    double A4 = Calc_A4(Pi4, t4);
                    double A5 = Calc_A5(Pi5, t5);
                    double A6 = Calc_A6(Pi6, t6);
                    double A7 = Calc_A7(Pi7, t7);
                    double A8 = Calc_A8(Pi8, t8);
                    double A9 = Calc_A9(Pi9, t9);
                    double A10 = Calc_A10(Pi10, t10);
                    double A11 = Calc_A11(Pi11, t11);
                    double A12 = Calc_A12(Pi12, t12);
                    double A13 = Calc_A13(Pi13, t13);
                    double ZnOv = Calc_ZnOv(A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13);

                    //==== Main equations for determination of steady state probabilities ========
                    P0[i, j] = Calc_P0(A0, ZnOv);
                    P1[i, j] = Calc_P1(A1, ZnOv);
                    P2[i, j] = Calc_P2(A2, ZnOv);
                    P3[i, j] = Calc_P3(A3, ZnOv);
                    P4[i, j] = Calc_P4(A4, ZnOv);
                    P5[i, j] = Calc_P5(A5, ZnOv);
                    P6[i, j] = Calc_P6(A6, ZnOv);
                    P7[i, j] = Calc_P7(A7, ZnOv);
                    P8[i, j] = Calc_P8(A8, ZnOv);
                    P9[i, j] = Calc_P9(A9, ZnOv);
                    P10[i, j] = Calc_P10(A10, ZnOv);
                    P11[i, j] = Calc_P11(A11, ZnOv);
                    P12[i, j] = Calc_P12(A12, ZnOv);
                    P13[i, j] = Calc_P13(A13, ZnOv);

                    AVL[i, j] = Calc_Avl(P0[i, j], P1[i, j], P2[i, j], P3[i, j], P4[i, j], P5[i, j], P6[i, j], P9[i, j], P11[i, j], P13[i, j]);
                    UNAVL[i, j] = Calc_Unavl(AVL[i, j]);
                }
            }
        }
        
        private int Calc_n(double Tmin, double Tmax, double Tstep)
        {
            int retVal = 0;

            do
            {
                retVal++;
                Tmin += Tstep;
            } while (Tmin <= Tmax && retVal <= n_max);

            return retVal;
        }

        private int Calc_m(double gamma_min, double gamma_max, double gamma_step)
        {
            decimal gamma_min_Dec = Convert.ToDecimal(gamma_min);
            decimal gamma_max_Dec = Convert.ToDecimal(gamma_max);
            decimal gamma_step_Dec = Convert.ToDecimal(gamma_step);

            int retVal = 0;

            do
            {
                retVal++;
                gamma_min_Dec += gamma_step_Dec;
            } while (gamma_min_Dec <= gamma_max_Dec && retVal <= m_max);

            return retVal;
        }

        private double Calc_mu1(double mu)
        {
            return mu;
        }

        private double Calc_mu2(double mu)
        {
            return mu / 2;
        }

        private double Calc_gamma(double gamma_max, double gamma_step, int iteration)
        {
            return (double)((decimal)gamma_max - (decimal)gamma_step * iteration);
        }

        private double Calc_gamma1(double k, double gamma)
        {
            return k * gamma;
        }

        private double Calc_gamma2(double gamma1)
        {
            return gamma1 / 2;
        }

        private double Calc_gamma3(double gamma2)
        {
            return gamma2;
        }

        private double Calc_T(double t_min, double t_step, int iteration)
        {
            return t_min + t_step * iteration;
        }

        private double Calc_t0(double gamma2, double T)
        {
            return 1 / gamma2 * (1 - Math.Exp(-gamma2 * T));
        }

        private double Calc_t1(double gamma1, double tcr)
        {
            return 1 / (2 * gamma1) * (1 - Math.Exp(-2 * gamma1 * tcr));
        }

        private double Calc_t2(double gamma2, double mu)
        {
            return Math.Pow(gamma2 / (gamma2 + mu), 2);
        }

        private double Calc_t3(double t2)
        {
            return t2;
        }

        private double Calc_t4(double gamma1, double tcr)
        {
            return 1 / gamma1 * (1 - Math.Exp(-gamma1 * tcr));
        }

        private double Calc_t5(double mu2)
        {
            return 5 / (4 * mu2);
        }

        private double Calc_t6(double t4)
        {
            return t4;
        }

        private double Calc_t7(double mu1)
        {
            return 2 / mu1;
        }

        private double Calc_t08(double T, double gamma2)
        {
            return 1 - (T * Math.Exp(-gamma2 * T) / (1 - gamma2 * T));
        }

        private double Calc_t8(double T, double t08)
        {
            return T - t08;
        }

        private double Calc_t9(double gamma3, double tcr)
        {
            return 1 / (2 * gamma3) * (1 - Math.Exp(-2 * gamma3 * tcr));
        }

        private double Calc_t10(double t8)
        {
            return t8;
        }

        private double Calc_t11(double gamma3, double tcr)
        {
            return 1 / gamma3 * (1 - Math.Exp(-gamma3 * tcr));
        }

        private double Calc_t12(double t8)
        {
            return t8;
        }

        private double Calc_t13(double t11)
        {
            return t11;
        }

        private double Calc_p01(double gamma2, double T)
        {
            return Math.Exp(-gamma2 * T);
        }

        private double Calc_p08(double p01)
        {
            return 1 - p01;
        }

        private double Calc_p10(double gamma1, double tcr)
        {
            return Math.Exp(-2 * gamma1 * tcr);
        }

        private double Calc_p12(double gamma1, double T)
        {
            return 1 / 2 * (1 - Math.Exp(-2 * gamma1 * T));
        }

        private double Calc_p13(double p10, double p12)
        {
            return 1 - p10 - p12;
        }

        private double Calc_p92(double gamma3, double T)
        {
            return 1 / 2 * (1 - Math.Exp(-2 * gamma3) * T);
        }

        private double Calc_p98(double gamma3, double tcr)
        {
            return Math.Exp(-2 * gamma3 * tcr);
        }

        private double Calc_p93(double p92, double p98)
        {
            return 1 - p92 - p98;
        }

        private double Calc_p24(double gamma2, double mu1, double T)
        {
            return (1 + mu1 * T) * Math.Exp(-(gamma2 + mu1) * T);
        }

        private double Calc_p34(double p24)
        {
            return p24;
        }

        private double Calc_p210(double gamma2, double mu1, double T)
        {
            return gamma2 * (1 + mu1 * T) * Math.Exp(-(gamma2 + mu1) * T);
        }

        private double Calc_p312(double p210)
        {
            return p210;
        }

        private double Calc_p20(double p24, double p210)
        {
            return 1 - p24 - p210;
        }

        private double Calc_p30(double p20)
        {
            return p20;
        }

        private double Calc_p42(double gamma1, double tcr)
        {
            return Math.Exp(-gamma1 * tcr);
        }

        private double Calc_p43(double p42)
        {
            return p42;
        }

        private double Calc_p65(double p42)
        {
            return p42;
        }

        private double Calc_p45(double p42, double p43)
        {
            return 1 - p42 - p43;
        }

        private double Calc_E1(double mu2, double T)
        {
            return Math.Exp(-2 * mu2 * T);
        }

        private double Calc_E2(double mu2, double T, double E1)
        {
            return 1 - Math.Pow(mu2, 2) * Math.Pow(T, 2) * E1 - E1;
        }

        private double Calc_p52(double mu2, double T, double E1, double E2)
        {
            return 1 / 2 * E2 - mu2 * T * E1;
        }

        private double Calc_p56(double mu2, double T)
        {
            return Math.Pow((1 + mu2 * T) * Math.Exp(-mu2 * T), 2);
        }

        private double Calc_p53(double p52, double p56)
        {
            return 1 - p52 - p56;
        }

        private double Calc_p67(double p65)
        {
            return 1 - p65;
        }

        private double Calc_p1110(double gamma3, double tcr)
        {
            return Math.Exp(-gamma3 * tcr);
        }

        private double Calc_p115(double p1110)
        {
            return 1 - p1110;
        }

        private double Calc_p1312(double p1110)
        {
            return p1110;
        }

        private double Calc_p135(double p115)
        {
            return p115;
        }

        private double Calc_a1(double p10)
        {
            return p10;
        }

        private double Calc_a2(double p20)
        {
            return p20;
        }

        private double Calc_a3(double p30)
        {
            return p30;
        }

        private double Calc_a4(double p01)
        {
            return p01;
        }

        private double Calc_a5(double p12)
        {
            return p12;
        }

        private double Calc_a6(double p42)
        {
            return p42;
        }

        private double Calc_a7(double p52)
        {
            return p52;
        }

        private double Calc_a8(double p92)
        {
            return p92;
        }

        private double Calc_a9(double p13)
        {
            return p13;
        }

        private double Calc_a10(double p43)
        {
            return p43;
        }

        private double Calc_a11(double p53)
        {
            return p53;
        }

        private double Calc_a12(double p93)
        {
            return p93;
        }

        private double Calc_a13(double p24)
        {
            return p24;
        }

        private double Calc_a14(double p34)
        {
            return p34;
        }

        private double Calc_a15(double p45)
        {
            return p45;
        }

        private double Calc_a16(double p65)
        {
            return p65;
        }

        private double Calc_a17(double p75)
        {
            return p75;
        }

        private double Calc_a18(double p115)
        {
            return p115;
        }

        private double Calc_a19(double p135)
        {
            return p135;
        }

        private double Calc_a20(double p56)
        {
            return p56;
        }

        private double Calc_a21(double p67)
        {
            return p67;
        }

        private double Calc_a22(double p08)
        {
            return p08;
        }

        private double Calc_a23(double p98)
        {
            return p98;
        }

        private double Calc_a24(double p210)
        {
            return p210;
        }

        private double Calc_a25(double p1110)
        {
            return p1110;
        }

        private double Calc_a26(double p312)
        {
            return p312;
        }

        private double Calc_a27(double p1312)
        {
            return p1312;
        }

        private double Calc_k7(double a20, double a21) { return a20 * a21; }

        private double Calc_k8(double a22, double a23)
        {
            return a22 / (1 - a23);
        }

        private double Calc_k10(double a24, double a25)
        {
            return a24 / (1 - a25);
        }

        private double Calc_k12(double a26, double a27)
        {
            return a26 / (1 - a27);
        }

        private double Calc_k20(double a4, double a5)
        {
            return a4 * a5;
        }

        private double Calc_k21(double a6, double a13)
        {
            return a6 * a13;
        }

        private double Calc_k22(double a6, double a14)
        {
            return a6 * a14;
        }

        private double Calc_k23(double a7, double a15)
        {
            return a7 * a15;
        }

        private double Calc_k28(double a8, double k8)
        {
            return a8 * k8;
        }

        private double Calc_k24(double a7, double a16, double a20)
        {
            return a7 * a16 * a20;
        }

        private double Calc_k25(double a7, double a17, double k7)
        {
            return a7 * a17 * k7;
        }

        private double Calc_k26(double a7, double a18, double k10)
        {
            return a7 * a18 * k10;
        }

        private double Calc_k27(double a7, double a19, double k12)
        {
            return a7 * a19 * k12;
        }

        private double Calc_k29(double a13, double k23)
        {
            return a13 * k23;
        }

        private double Calc_k30(double a14, double k23)
        {
            return a14 * k23;
        }

        private double Calc_k31(double a4, double a9)
        {
            return a4 * a9;
        }

        private double Calc_k32(double a10, double a13)
        {
            return a10 * a13;
        }

        private double Calc_k33(double a10, double a14)
        {
            return a10 * a14;
        }

        private double Calc_k34(double a11, double a15)
        {
            return a11 * a15;
        }

        private double Calc_k35(double a11, double a16, double a20)
        {
            return a11 * a16 * a20;
        }

        private double Calc_k36(double a11, double a17, double k7)
        {
            return a11 * a17 * k7;
        }

        private double Calc_k37(double a11, double a18, double k10)
        {
            return a11 * a18 * k10;
        }

        private double Calc_k38(double a11, double a19, double k12)
        {
            return a11 * a19 * k12;
        }

        private double Calc_k39(double a12, double k8)
        {
            return a12 * k8;
        }

        private double Calc_k50(double a16, double a20)
        {
            return a16 * a20;
        }

        private double Calc_k51(double a17, double k7)
        {
            return a17 * k7;
        }

        private double Calc_k52(double a18, double k10)
        {
            return a18 * k10;
        }

        private double Calc_k53(double a19, double k12)
        {
            return a19 * k12;
        }

        private double Calc_k54(double a13, double a15)
        {
            return a13 * a15;
        }

        private double Calc_k55(double a14, double a15)
        {
            return a14 * a15;
        }

        private double Calc_Zn5(double k50, double k51)
        {
            return 1 - k50 - k51;
        }

        private double Calc_c20(double k24, double k52)
        {
            return k24 * k52;
        }

        private double Calc_c21(double k24, double k54)
        {
            return k24 * k54;
        }

        private double Calc_c22(double k24, double k53)
        {
            return k24 * k53;
        }

        private double Calc_c23(double k24, double k55)
        {
            return k24 * k55;
        }

        private double Calc_c24(double k25, double k52)
        {
            return k25 * k52;
        }

        private double Calc_c25(double k25, double k54)
        {
            return k25 * k54;
        }

        private double Calc_c26(double k25, double k53)
        {
            return k25 * k53;
        }

        private double Calc_c27(double k25, double k55)
        {
            return k25 * k55;
        }

        private double Calc_c28(double k20, double Zn5)
        {
            return k20 * Zn5;
        }

        private double Calc_c29(double k28, double Zn5)
        {
            return k28 * Zn5;
        }

        private double Calc_c30(double k22, double Zn5)
        {
            return k22 * Zn5;
        }

        private double Calc_c31(double k27, double Zn5)
        {
            return k27 * Zn5;
        }

        private double Calc_c32(double k30, double Zn5)
        {
            return k30 * Zn5;
        }

        private double Calc_c33(double a13, double c28)
        {
            return a13 * c28;
        }

        private double Calc_c34(double a13, double c29)
        {
            return a13 * c29;
        }

        private double Calc_c35(double a13, double c22)
        {
            return a13 * c22;
        }

        private double Calc_c36(double a13, double c23)
        {
            return a13 * c23;
        }

        private double Calc_c37(double a13, double c26)
        {
            return a13 * c26;
        }

        private double Calc_c38(double a13, double c27)
        {
            return a13 * c27;
        }

        private double Calc_c39(double a13, double c30)
        {
            return a13 * c30;
        }

        private double Calc_c40(double a13, double c31)
        {
            return a13 * c31;
        }

        private double Calc_c41(double a13, double c32)
        {
            return a13 * c32;
        }

        private double Calc_Zn2(double Zn5, double c20, double c21, double c24, double c25, double k21, double k26, double k29)
        {
            return Zn5 * (1 - c20 - c21 - c24 - c25 - k21 - k26 - k29);
        }

        private double Calc_c42(double a14, double Zn2)
        {
            return a14 * Zn2;
        }

        private double Calc_c43(double k52, double c28)
        {
            return k52 * c28;
        }

        private double Calc_c44(double k52, double c29)
        {
            return k52 * c29;
        }

        private double Calc_c45(double k52, double c22)
        {
            return k52 * c22;
        }

        private double Calc_c46(double k52, double c23)
        {
            return k52 * c23;
        }

        private double Calc_c47(double k52, double c26)
        {
            return k52 * c26;
        }

        private double Calc_c48(double k52, double c27)
        {
            return k52 * c27;
        }

        private double Calc_c49(double k52, double c30)
        {
            return k52 * c30;
        }

        private double Calc_c50(double k52, double c31)
        {
            return k52 * c31;
        }

        private double Calc_c51(double k52, double c32)
        {
            return k52 * c32;
        }

        private double Calc_c52(double k54, double c28)
        {
            return k54 * c28;
        }

        private double Calc_c53(double k54, double c29)
        {
            return k54 * c29;
        }

        private double Calc_c54(double k54, double c22)
        {
            return k54 * c22;
        }

        private double Calc_c55(double k54, double c23)
        {
            return k54 * c23;
        }

        private double Calc_c56(double k54, double c26)
        {
            return k54 * c26;
        }

        private double Calc_c57(double k54, double c27)
        {
            return k54 * c27;
        }

        private double Calc_c58(double k54, double c30)
        {
            return k54 * c30;
        }

        private double Calc_c59(double k54, double c31)
        {
            return k54 * c31;
        }

        private double Calc_c60(double k54, double c32)
        {
            return k54 * c32;
        }

        private double Calc_c61(double k53, double Zn2)
        {
            return k53 * Zn2;
        }

        private double Calc_c62(double k55, double Zn2)
        {
            return k55 * Zn2;
        }

        private double Calc_kf41(double c33, double c34)
        {
            return c33 + c34;
        }

        private double Calc_kf42(double c35, double c36, double c37, double c38, double c39, double c40, double c41, double c42)
        {
            return c35 + c36 + c37 + c38 + c39 + c40 + c41 + c42;
        }

        private double Calc_kf51(double c43, double c44, double c52, double c53)
        {
            return c43 + c44 + c52 + c53;
        }

        private double Calc_kf52(double c45, double c46, double c47, double c48, double c49, double c50, double c51, double c54, double c55,
            double c56, double c57, double c58, double c59, double c60, double c61, double c62)
        {
            return c45 + c46 + c47 + c48 + c49 + c50 + c51 + c54 + c55 + c56 + c57 + c58 + c59 + c60 + c61 + c62;
        }

        private double Calc_kf300(double Zn2, double Zn5)
        {
            return Zn2 * Zn5;
        }

        private double Calc_kf301(double k31, double Zn2, double Zn5)
        {
            return k31 * Zn2 * Zn5;
        }

        private double Calc_kf302(double k39, double Zn2, double Zn5)
        {
            return k39 * Zn2 * Zn5;
        }

        private double Calc_kf303(double k32, double c28, double Zn5)
        {
            return k32 * c28 * Zn5;
        }

        private double Calc_kf304(double k32, double c29, double Zn5)
        {
            return k32 * c29 * Zn5;
        }

        private double Calc_kf305(double k32, double c22, double Zn5)
        {
            return k32 * c22 * Zn5;
        }

        private double Calc_kf306(double k32, double c23, double Zn5)
        {
            return k32 * c23 * Zn5;
        }

        private double Calc_kf307(double k32, double c26, double Zn5)
        {
            return k32 * c26 * Zn5;
        }

        private double Calc_kf308(double k32, double c27, double Zn5)
        {
            return k32 * c27 * Zn5;
        }

        private double Calc_kf309(double k32, double c30, double Zn5)
        {
            return k32 * c30 * Zn5;
        }

        private double Calc_kf310(double k32, double c31, double Zn5)
        {
            return k32 * c31 * Zn5;
        }

        private double Calc_kf311(double k32, double c32, double Zn5)
        {
            return k32 * c32 * Zn5;
        }

        private double Calc_kf312(double k37, double c28, double Zn5)
        {
            return k37 * c28 * Zn5;
        }

        private double Calc_kf313(double k37, double c29, double Zn5)
        {
            return k37 * c29 * Zn5;
        }

        private double Calc_kf314(double k37, double c22, double Zn5)
        {
            return k37 * c22 * Zn5;
        }

        private double Calc_kf315(double k37, double c23, double Zn5)
        {
            return k37 * c23 * Zn5;
        }

        private double Calc_kf316(double k37, double c26, double Zn5)
        {
            return k37 * c26 * Zn5;
        }

        private double Calc_kf317(double k37, double c27, double Zn5)
        {
            return k37 * c27 * Zn5;
        }

        private double Calc_kf318(double k37, double c30, double Zn5)
        {
            return k37 * c30 * Zn5;
        }

        private double Calc_kf319(double k37, double c31, double Zn5)
        {
            return k37 * c31 * Zn5;
        }

        private double Calc_kf320(double k37, double c32, double Zn5)
        {
            return k37 * c32 * Zn5;
        }

        private double Calc_kf321(double k33, double Zn2, double Zn5)
        {
            return k33 * Zn2 * Zn5;
        }

        private double Calc_kf322(double k38, double Zn2, double Zn5)
        {
            return k38 * Zn2 * Zn5;
        }

        private double Calc_kf323(double k34, double kf41, double Zn5)
        {
            return k34 * kf41 * Zn5;
        }

        private double Calc_kf324(double k34, double kf42, double Zn5)
        {
            return k34 * kf42 * Zn5;
        }

        private double Calc_kf325(double k35, double kf51)
        {
            return k35 * kf51;
        }

        private double Calc_kf326(double k35, double kf52)
        {
            return k35 * kf52;
        }

        private double Calc_kf327(double kf325)
        {
            return kf325;
        }

        private double Calc_kf328(double k36, double kf52)
        {
            return k36 * kf52;
        }

        private double Calc_Ch3(double kf301, double kf302, double kf303, double kf304, double kf312, double kf313,
            double kf323, double kf325, double kf327)
        {
            return kf301 + kf302 + kf303 + kf304 + kf312 + kf313 + kf323 + kf325 + kf327;
        }

        private double Calc_Zn3(double kf300, double kf305, double kf306, double kf307, double kf308, double kf309, double kf310, double kf311, double kf314,
            double kf315, double kf316, double kf317, double kf318, double kf319, double kf320, double kf321, double kf322, double kf324, double kf326, double kf328)
        {
            return kf300 - kf305 - kf306 - kf307 - kf308 - kf309 - kf310 - kf311 - kf314 - kf315 - kf316 -
                kf317 - kf318 - kf319 - kf320 - kf321 - kf322 - kf324 - kf326 - kf328;
        }

        private double Calc_kf200(double Zn2, double Zn3)
        {
            return Zn2 * Zn3;
        }

        private double Calc_kf201(double c28, double Zn3)
        {
            return c28 * Zn3;
        }

        private double Calc_kf202(double c29, double Zn3)
        {
            return c29 * Zn3;
        }

        private double Calc_kf203(double c22, double Ch3)
        {
            return c22 * Ch3;
        }

        private double Calc_kf204(double c23, double Ch3)
        {
            return c23 * Ch3;
        }

        private double Calc_kf205(double c26, double Ch3)
        {
            return c26 * Ch3;
        }

        private double Calc_kf206(double c27, double Ch3)
        {
            return c27 * Ch3;
        }

        private double Calc_kf207(double c30, double Ch3)
        {
            return c30 * Ch3;
        }

        private double Calc_kf208(double c31, double Ch3)
        {
            return c31 * Ch3;
        }

        private double Calc_kf209(double c32, double Ch3)
        {
            return c32 * Ch3;
        }

        private double Calc_Ch2(double kf201, double kf202, double kf203, double kf204, double kf205, double kf206, double kf207, double kf208, double kf209)
        {
            return kf201 + kf202 + kf203 + kf204 + kf205 + kf206 + kf207 + kf208 + kf209;
        }

        private double Calc_kf100(double k10, double Ch2, double kf200)
        {
            return (k10 * Ch2) / kf200;
        }

        private double Calc_kf120(double k12, double kf300)
        {
            return k12 * kf300;
        }

        private double Calc_kf400(double a13, double Ch2, double Zn3)
        {
            return a13 * Ch2 * Zn3;
        }

        private double Calc_kf401(double a14, double Ch3, double kf200)
        {
            return a14 * Ch3 * kf200;
        }

        private double Calc_kf402(double kf400, double kf401, double kf200, double Zn3)
        {
            return (kf400 + kf401) / (kf200 * Zn3);
        }

        private double Calc_kf500(double kf51, double Zn3)
        {
            return kf51 * Zn3;
        }

        private double Calc_kf501(double kf52, double Ch3)
        {
            return kf52 * Ch3;
        }

        private double Calc_kf502(double Zn2, double Zn3, double Zn5)
        {
            return Zn2 * Zn3 * Zn5;
        }

        private double Calc_kf503(double kf500, double kf501, double kf502)
        {
            return (kf500 + kf501) / kf502;
        }

        private double Calc_kf600(double a20, double kf503)
        {
            return a20 * kf503;
        }

        private double Calc_kf700(double a21, double kf600)
        {
            return a21 * kf600;
        }

        private double Calc_kf800(double a22, double k8)
        {
            return a22 * k8;
        }

        private double Calc_kf201_1(double Ch2, double kf200)
        {
            return Ch2 / kf200;
        }

        private double Calc_U(double a4, double kf201_1, double kf300, double kf402,
            double kf503, double kf600, double kf700, double kf800, double kf100, double kf120)
        {
            return 1 + a4 + kf201_1 + kf300 + kf402 + kf503 + kf600 + kf700 + 2 * (kf800 + kf100 + kf120);
        }

        private double Calc_Pi0(double U)
        {
            return 1 / U;
        }

        private double Calc_Pi1(double a4, double Pi0)
        {
            return a4 * Pi0;
        }

        private double Calc_Pi2(double kf201_1, double Pi0)
        {
            return kf201_1 * Pi0;
        }

        private double Calc_Pi3(double kf300, double Pi0)
        {
            return kf300 * Pi0;
        }

        private double Calc_Pi4(double kf402, double Pi0)
        {
            return kf402 * Pi0;
        }

        private double Calc_Pi5(double kf503, double Pi0)
        {
            return kf503 * Pi0;
        }

        private double Calc_Pi6(double kf600, double Pi0)
        {
            return kf600 * Pi0;
        }

        private double Calc_Pi7(double kf700, double Pi0)
        {
            return kf700 * Pi0;
        }

        private double Calc_Pi8(double kf800, double Pi0)
        {
            return kf800 * Pi0;
        }

        private double Calc_Pi9(double Pi8)
        {
            return Pi8;
        }

        private double Calc_Pi10(double kf100, double Pi0)
        {
            return kf100 * Pi0;
        }

        private double Calc_Pi11(double Pi10)
        {
            return Pi10;
        }

        private double Calc_Pi12(double kf120, double Pi0)
        {
            return kf120 * Pi0;
        }

        private double Calc_Pi13(double Pi12)
        {
            return Pi12;
        }

        private double Calc_A0(double Pi0, double t0)
        {
            return Pi0 * t0;
        }

        private double Calc_A1(double Pi1, double t1)
        {
            return Pi1 * t1;
        }

        private double Calc_A2(double Pi2, double t2)
        {
            return Pi2 * t2;
        }

        private double Calc_A3(double Pi3, double t3)
        {
            return Pi3 * t3;
        }

        private double Calc_A4(double Pi4, double t4)
        {
            return Pi4 * t4;
        }

        private double Calc_A5(double Pi5, double t5)
        {
            return Pi5 * t5;
        }

        private double Calc_A6(double Pi6, double t6)
        {
            return Pi6 * t6;
        }

        private double Calc_A7(double Pi7, double t7)
        {
            return Pi7 * t7;
        }

        private double Calc_A8(double Pi8, double t8)
        {
            return Pi8 * t8;
        }

        private double Calc_A9(double Pi9, double t9)
        {
            return Pi9 * t9;
        }

        private double Calc_A10(double Pi10, double t10)
        {
            return Pi10 * t10;
        }

        private double Calc_A11(double Pi11, double t11)
        {
            return Pi11 * t11;
        }

        private double Calc_A12(double Pi12, double t12)
        {
            return Pi12 * t12;
        }

        private double Calc_A13(double Pi13, double t13)
        {
            return Pi13 * t13;
        }

        private double Calc_ZnOv(double A0, double A1, double A2, double A3, double A4, double A5, double A6,
            double A7, double A8, double A9, double A10, double A11, double A12, double A13)
        {
            return A0 + A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9 + A10 + A11 + A12 + A13;
        }

        private double Calc_P0(double A0, double ZnOv)
        {
            return A0 / ZnOv;
        }

        private double Calc_P1(double A1, double ZnOv)
        {
            return A1 / ZnOv;
        }

        private double Calc_P2(double A2, double ZnOv)
        {
            return A2 / ZnOv;
        }

        private double Calc_P3(double A3, double ZnOv)
        {
            return A3 / ZnOv;
        }

        private double Calc_P4(double A4, double ZnOv)
        {
            return A4 / ZnOv;
        }

        private double Calc_P5(double A5, double ZnOv)
        {
            return A5 / ZnOv;
        }

        private double Calc_P6(double A6, double ZnOv)
        {
            return A6 / ZnOv;
        }

        private double Calc_P7(double A7, double ZnOv)
        {
            return A7 / ZnOv;
        }

        private double Calc_P8(double A8, double ZnOv)
        {
            return A8 / ZnOv;
        }

        private double Calc_P9(double A9, double ZnOv)
        {
            return A9 / ZnOv;
        }

        private double Calc_P10(double A10, double ZnOv)
        {
            return A10 / ZnOv;
        }

        private double Calc_P11(double A11, double ZnOv)
        {
            return A11 / ZnOv;
        }

        private double Calc_P12(double A12, double ZnOv)
        {
            return A12 / ZnOv;
        }

        private double Calc_P13(double A13, double ZnOv)
        {
            return A13 / ZnOv;
        }

        private double Calc_Avl(double P0, double P1, double P2, double P3, double P4, double P5, double P6, double P9, double P11, double P13)
        {
            return P0 + P1 + P2 + P3 + P4 + P5 + P6 + P9 + P11 + P13;
        }

        private double Calc_Unavl(double avl)
        {
            return 1 - avl;
        }
    }
}