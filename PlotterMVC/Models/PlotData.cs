using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace PlotterMVC.Models
{
    public class PlotData
    {
        public double[] T { get; set; }
        public double[] gamma { get; set;}
        public double[,] A { get; set; }
    }
}

