using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Web;

namespace PlotterMVC.Globalization
{
    public class NumberProvider
    {
        public static NumberFormatInfo GetFormatter()
        {
            NumberFormatInfo nInfo = new NumberFormatInfo();
            nInfo.NumberDecimalSeparator = ".";

            return nInfo;
        }
    }
}