using System;
using System.IO;
using System.Collections.Generic;

namespace Fulleren2
{
    class Program
    {
        static double d = 1.42;
        static double a = 100;
        static double eps = 0.5;
        static double f = 0.5;
        static double b = 1;
        static double dx = 0.01;
        static int MAXR = 1000;//расстояние между атомами, значительно превышающее остальные расстояния

        class Atom
        {
            public double x, y, z;
            public char type;
            public Atom(double x, double y, double z, char type)
            {
                this.x = x;
                this.y = y;
                this.z = z;
                this.type = type;
            }
            public Atom(Atom atom)
            {
                this.x = atom.x;
                this.y = atom.y;
                this.z = atom.z;
                this.type = atom.type;
            }
        };
        class Graphene
        {
            public Atom[] atoms;
            public int N;
            public Graphene(List<Atom> list)
            {
                N = list.Count;
                atoms = new Atom[list.Count];
                int i = 0;
                foreach (Atom atom in list)
                {
                    atoms[i] = new Atom(atom);
                    i++;
                }
            }
            public Graphene(Graphene graphene, int i, double dx, double dy, double dz)
            {
                this.N = graphene.N;
                this.atoms = graphene.atoms;
                this.atoms[i].x += dx;
                this.atoms[i].y += dy;
                this.atoms[i].z += dz;
            }
            public double R(int i, int j)
            {
                if (i != j)
                    return (Math.Sqrt(Math.Pow((atoms[i].x - atoms[j].x), 2) + Math.Pow((atoms[i].y - atoms[j].y), 2) + Math.Pow((atoms[i].z - atoms[j].z), 2)));
                else return (MAXR);
            }
            public double ANG(int i, int j, int k)
            {
                return (180.0 / Math.PI * Math.Acos((R(i, j) * R(i, j) + R(j, k) * R(j, k) - R(i, k) * R(i, k)) / (2 * R(i, j) * R(j, k))));
            }
            public double Energy1()
            {
                double U = 0;
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        if (Math.Abs(R(i, j) - d) < eps)
                        {
                            U = U + a * Math.Pow(Math.Abs(R(i, j) - d), 2);
                        }
                    }
                }
                return U;
            }
            public double Energy2(int[,] mat)
            {
                double U = 0;
                for (int i = 0; i < N; i++)
                {
                    U = U + b * Math.Pow(Math.Abs((ANG(mat[i, 0] - 1, i, mat[i, 1] - 1)) - 120), 2);

                    if (mat[i, 2] != 0)
                    {
                        U = U + b * Math.Pow(Math.Abs((ANG(mat[i, 0] - 1, i, mat[i, 2] - 1)) - 120), 2);
                        U = U + b * Math.Pow(Math.Abs((ANG(mat[i, 1] - 1, i, mat[i, 2] - 1)) - 120), 2);
                    }
                }
                return U;
            }
        };

        static double OLS(double[] f, double[] E)
        {
            int N = f.Length;
            double[] f2 = new double[f.Length];
            for (int i = 0; i < N; i++)
                f2[i] = f[i] * f[i];
            double m11, m12, m21, m22, g1, g2;
            m11 = N;
            m12 = m21 = m22 = g1 = g2 = 0;
            for (int i = 0; i < N; i++)
            {
                m12 += f2[i];
                m22 += f2[i] * f2[i];
                g1 += E[i];
                g2 += E[i] * f2[i];
            }
            m21 = m12;
            return (g1 - m11 / m12 * g2) / (m12 - m11 / m21 * m22);
                
        }

        static void Main(string[] args)
        {


            List<Atom> list = new List<Atom>();

            using (StreamReader sr = new StreamReader("input.txt"))
            {
                string numbers;
                while ((numbers = sr.ReadLine()) != null)
                {
                    string[] arr = numbers.Split('\t');
                    double x = Double.Parse(arr[0]);
                    double y = Double.Parse(arr[1]);
                    double z = Double.Parse(arr[2]);
                    char type = Char.Parse(arr[3]);
                    Console.WriteLine(x + "\t" + y + "\t" + z + "\t" + type);
                    list.Add(new Atom(x, y, z, type));
                }
            }

            Graphene graphene = new Graphene(list);
            int N = list.Count;
            /*Заполнение матрицы соседей*/
            int[,] mat = new int[N, 3];
            double[] r = new double[N];
            for(int i = 0; i < N; i++)
            {
                int count = 0;
                for(int j = 0; j < N; j++)
                {
                    if (Math.Abs(graphene.R(i, j) - d) < f)
                    {
                        mat[i, count] = j + 1;
                        count++;
                    }
                }
            }
            Console.WriteLine();
            for (int i = 0; i < N; i++)
            {
                Console.Write(i + 1 + "\t");
                for (int j = 0; j < 3; j++)
                    Console.Write(mat[i, j] + "\t");
                Console.WriteLine();
            }

            /*Запись файла*/
            using (StreamWriter sw = new StreamWriter("output.txt"))
            {
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < 3; j++)
                        sw.Write(mat[i, j] + " ");
            }

            /*Подсчет потенциала*/
            double E = 0;
            E += graphene.Energy1();
            Console.WriteLine("E = " + E);
            E += graphene.Energy2(mat);
            Console.WriteLine("E = " + E);
            /*
            Graphene graphene1 = new Graphene(graphene, 0, dx, 0, 0);
            double E1 = 0;
            E1 += graphene1.Energy1();
            Console.WriteLine("E1 = " + E1);
            E1 += graphene1.Energy2(mat);
            Console.WriteLine("E1 = " + E1);
            double f1 = (E - E1) / dx;
            Console.WriteLine("f1 = " + f1);*/
            double[] f1 = new double[] { -2, -1, 0, 1, 3, 4 };
            double[] E1 = new double[] { 3, 1, 1, 0, -1, -3 };
            Console.WriteLine("k = " + OLS(f1, E1)); 

            Console.ReadLine();

        }
    }
}
