using PlotterMVC.Models;
using PlotterMVC.SMM;
using System;
using System.Reflection;
using System.Collections.Generic;
using System.Linq;
using System.Web;
using System.Web.Mvc;
using System.Web.Script.Serialization;
using System.Globalization;
using PlotterMVC.Globalization;

namespace PlotterMVC.Controllers
{
    public class HomeController : Controller
    {
        // GET: Home
        public ActionResult Index()
        {
            PlotInitials initData = GetInitialData();

            ViewBag.Tmin = initData.Tmin;
            ViewBag.Tmax = initData.Tmax;
            ViewBag.mu = initData.mu;
            ViewBag.gamma_min = Convert.ToDecimal(initData.gamma_min).ToString(NumberProvider.GetFormatter());
            ViewBag.gamma_max = Convert.ToDecimal(initData.gamma_max).ToString(NumberProvider.GetFormatter());

            return View();
        }
        
        private PlotInitials GetInitialData()
        {
            PlotInitials initData = new PlotInitials();

            IaaS iaas = new IaaS();
            
            initData.Tmin = iaas.t_min;
            initData.Tmax = iaas.t_max;
            initData.gamma_min = iaas.gamma_min;
            initData.gamma_max = iaas.gamma_max;
            initData.mu = iaas.mu;

            return initData;
        }

        [HttpPost]
        public string GetData(string args)
        {
            JavaScriptSerializer jsSerializer = new JavaScriptSerializer();
            PlotInitials argsAsObj = jsSerializer.Deserialize<PlotInitials>(args);

            IaaS iaas = new IaaS();

            iaas.Set_T_min(argsAsObj.Tmin);
            iaas.Set_T_max(argsAsObj.Tmax);
            iaas.Set_gamma_max(argsAsObj.gamma_max);
            iaas.Set_gamma_min(argsAsObj.gamma_min);
            iaas.Set_mu(argsAsObj.mu);

            iaas.Calculate();

            PlotData plotData = new PlotData();
            plotData.T = iaas.T;
            plotData.gamma = iaas.gamma;
            plotData.A = iaas.AVL;

            return jsSerializer.Serialize(plotData);
        }
    }
}