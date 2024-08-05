#### Please do not run this file directly as it will not work.
#### Please refer to file 00_run.R

# require(odeintr)

mycompile_sys <- function (name, sys, pars = NULL, const = FALSE, method = "rk5_i",
    sys_dim = -1L, atol = 1e-06, rtol = 1e-06, globals = "",
    headers = "", footers = "", compile = TRUE, observer = NULL,
    env = new.env(), ...){

    if (!is.null(pars)) {
        pcode = JLhandle_pars(pars, const)
        globals = paste0(pcode$globals, globals, collapse = ";\n")
        if (!is.null(pcode$getter))
            footers = paste0(pcode$getter, footers, collapse = "\n")
        if (!is.null(pcode$setter))
            footers = paste0(pcode$setter, footers, collapse = "\n")

    }
    sys = paste0(sys, collapse = "; ")
    stepper = make_stepper_type(method)
    stepper_constr = make_stepper_constr(method, atol, rtol)
    if (ceiling(sys_dim) < 1)
        sys_dim = get_sys_dim(sys)
    if (is.na(sys_dim)) {
        sys = vectorize_1d_sys(sys)
        sys_dim = 1L
    }
    code = if (!is.null(observer))
        read_template("compile_sys_r_obs_template")
    else read_template("compile_sys_template")
    code = gsub("__STEPPER_TYPE__", stepper, code)
    code = gsub("__STEPPER_CONSTRUCT__", stepper_constr, code)
    code = sub("__SYS_SIZE__", ceiling(sys_dim), code)
    code = sub("__GLOBALS__", globals, code)
    code = sub("__SYS__", sys, code)
    code = sub("__HEADERS__", headers, code)
    code = sub("__FOOTERS__", footers, code)
    code = gsub("__FUNCNAME__", name, code)

    sink(file("compiling_output.txt"), append = FALSE, type = "output") # Writing console output to log file

    if (compile) {
        Sys.setenv(CXX_STD = "CXX11")

        # sink(file("compiling_output.txt"), append = FALSE, type = "output") # Writing console output to log file
        # res = try(Rcpp::sourceCpp(code = code, env = env, ...))
        # closeAllConnections() # Close connection to log file      

        res = try(Rcpp::sourceCpp(code = code, verbose=TRUE, showOutput=TRUE, rebuild=TRUE, echo=FALSE, env = env, ...))
    
        if (!inherits(res, "try-error")) {
            if (!is.null(observer)) {
                do.call(paste0(name, "_set_observer"), list(f = observer),
                envir = env)
              do.call(paste0(name, "_set_output_processor"),
                list(f = proc_output), envir = env)
          }
          if (!identical(env, globalenv())) {
              if (name %in% search())
                detach(pos = match(name, search()))
              do.call("attach", list(what = env, name = name))
          }
      }
  }
          closeAllConnections() # Close connection to log file      

  return(invisible(code))
}

read_template = function(name)
{
  tn = paste0(name, ".cpp")
  fn = system.file(file.path("templates", tn),
                   package = "odeintr", mustWork = TRUE)
  con = file(fn)
  res = readLines(con)
  close(con)
  paste0(res, collapse = "\n")
}


# vectorize_1d_sys = function(sys)
# {
#   if (grepl("\\[\\.*\\]", sys)) return(sys)
#   sys = gsub("\\bx\\b", "x[0]", sys)
#   sys = gsub("\\bdxdt\\b", "dxdt[0]", sys)
#   return(sys)
# }

make_stepper_constr = function(method, atol, rtol)
{
  stepper_constr = "stepper_type()"
  if (grepl("euler_|rk4_|rk54_i|rk78_i|bs_|bsd_", method))
    stop("Invalid integration method")
  if (grepl("_a$", method))
    stepper_constr = paste0("odeint::make_controlled(", atol, ", ", rtol, ", stepper_type())")
  if (grepl("_i$", method))
    stepper_constr = paste0("odeint::make_dense_output(", atol, ", ", rtol, ", stepper_type())")
  return(stepper_constr)
}

# make_stepper_type_openmp = function(stepper)
# {
#   stepper = make_stepper_type(stepper)
#   sub(">", ", double, state_type, double, odeint::openmp_range_algebra>", stepper)
# }

make_stepper_type = function(stepper)
{
  stepper = sub("_i$|_a$", "", stepper)
  if (grepl("ab|am|abm", stepper))
  {
    steps = as.integer(sub("ab|am|abm([0-9]+)", "\\1", stepper))
    stepper = sub("(ab|am|abm)[0-9]+", "\\1", stepper)
  }
  switch(stepper,
         euler = "euler<state_type>",
         rk4 = "runge_kutta4<state_type>",
         rk54 = "runge_kutta_cash_karp54<state_type>",
         rk5 = "runge_kutta_dopri5<state_type>",
         rk78 = "runge_kutta_fehlberg78<state_type>",
         ab = paste0("adams_bashforth<", steps, ", state_type>"),
         am = paste0("adams_moulton<", steps, ", state_type>"),
         abm = paste0("adams_bashforth_moulton<", steps, ", state_type>"),
         bs = "bulirsch_stoer<state_type>",
         bsd = "bulirsch_stoer_dense_out<state_type>",
         paste0(stepper, "<state_type>"))
}

# get_sys_dim = function(x)
# {
#   x = gsub("\\[\\s*(\\d+)\\s*\\]", "\\[\\1\\]", x)
#   matches = gregexpr("dxdt\\[\\d+\\]", x)
#   lens = attr(matches[[1]], "match.length") - 7L
#   starts = unlist(matches) + 5L
#   indices = rep(NA, length(starts))
#   for (i in seq(along = indices))
#     indices[i] = substr(x, starts[i], starts[i] + lens)
#   return(max(as.integer(indices)) + 1L)
# }

# disable_asserts = function(makevars)
# {
#   con = pipe(paste("R CMD config CPPFLAGS"))
#   flags = readLines(con); close(con)
#   if (!is.finite(pmatch("-DNDEBUG", flags)))
#     flags = paste(flags, "-DNDEBUG")
#   flags = gsub("^\\s+|\\s+$", "", flags)
#   cat(paste0("CPPFLAGS=", flags, "\n"),
#       file = makevars, append = TRUE)
# }

# substitute_opt_level = function(flags, level, omit.debug)
# {
#   flags = gsub("-O\\d+", paste0("-O", level), flags)
#   if (omit.debug) flags = gsub("\\s*-g\\s*", "", flags)
#   flags = gsub("^\\s+|\\s+$", "", flags)
#   return(flags)
# }

# process_flags = function(name, level, omit.debug)
# {
#   con = pipe(paste("R CMD config", name))
#   flags = readLines(con); close(con)
#   flags = substitute_opt_level(flags, level, omit.debug)
#   paste0(name, "=", flags)
# }

# Jacobian1 = function(f)
# {
#   sep = ".."
#   vn = names(formals(f))[1]
#   e = body(f)
#   for (i in length(e))
#   {
#     es = deparse(e[[i]])
#     if (grepl(paste0("\\b", vn, "\\b"), es))
#     {
#       es = gsub("\\[\\s*(\\d+)\\s*\\]", paste0(sep, "\\1"), es)
#       de = deparse(stats::D(parse(text = es), vn))
#       de = gsub(paste0("\\bx", sep, "(\\d+)"), "x\\[\\1\\]", de)
#       e[[i]] = parse(text = de)
#     }
#   }
#   res = function() NULL
#   formals(res) = formals(f)
#   body(res) = e
#   return(res)
# }

# Jacobian2 = function(code, sys_dim = -1)
# {
#   sep = ".."
#   code = paste0(code, collapse = ";")
#   if (sys_dim < 1)
#     sys_dim = get_sys_dim(code)
#   if (is.na(sys_dim))
#   {
#     code = vectorize_1d_sys(code)
#     sys_dim = 1L
#   }
#   code = gsub("^\\s+|\\s+$", "", unlist(strsplit(code, ";")))
#   code = code[nzchar(code) != 0]
#   i = unlist(lapply(code, function(x)
#     as.numeric(sub("\\bdxdt\\[\\s*(\\d+)\\s*\\].*", "\\1", x))))
#   code = code[order(i)]
#   g = function(j, i, rhs)
#   {
#     var = paste0("x", sep, j - 1)
#     deriv = stats::D(parse(text = rhs), var)
#     deriv = deparse(deriv)
#     deriv = gsub(paste0("\\bx", sep, "(\\d+)"), "x\\[\\1\\]", deriv)
#     deriv = paste0("J(", i - 1, ", ", j - 1, ") = ", deriv, ";")
#   }
#   f = function(i)
#   {
#     rhs = sub("\\bdxdt\\[\\s*\\d+\\s*\\]\\s*=\\s*(.*)", "\\1", code[i])
#     rhs = gsub("\\[\\s*(\\d+)\\s*\\]", paste0(sep, "\\1"), rhs)
#     return(unlist(lapply(1:sys_dim, g, i, rhs)))
#   }
#   paste0(unlist(lapply(1:sys_dim, f)), collapse = "\n")
# }

###JL
JLhandle_pars = function(pars,  const = FALSE)
{
  switch(mode(pars),
         numeric =
         {
           stop("JL: ive changed the original code to not allow numeric input. See code.")
         },
         character =
           list(globals = JLmake_pars_decl(pars),
                setter = JLmake_pars_setter(pars),
                getter = JLmake_pars_getter(pars)),
         stop("Invalid parameter specification"))
}


# make_init_pars_decl = function(pars, const = FALSE)
# {
#   res = paste(names(pars), "=", pars, collapse = ", ")
#   res = paste0("double ", res, ";")
#   if (const) res = paste("const", res)
#   res
# }

# make_init_pars_setter = function(pars)
# {
#   pn = names(pars)
#   body = paste0("odeintr::", pn, " = ", pn, collapse = ";\n")
#   pars = paste0("double ", pn, " = ", pars, collapse = ", ")
#   code = read_template("pars_setter_template")
#   code = sub("__PARS__", pars, code)
#   code = sub("__BODY__", body, code)
#   return(paste0(code, collapse = "\n"))
# }

# make_init_pars_getter = function(pars)
# {
#   pn = names(pars)
#   body = paste0("out[\"", pn, "\"] = odeintr::", pn, collapse = ";\n")
#   code = read_template("pars_getter_template")
#   code = sub("__BODY__", body, code)
#   return(paste0(code, collapse = "\n"))
# }

# make_vec_pars_decl = function(N)
# {
#   setter = read_template("pars_vec_setter_template")
#   getter = read_template("pars_vec_getter_template")
#   list(globals = paste0("std::array<double, ", N, "> pars;"),
#        setter = paste0(setter, collapse = "\n"),
#        getter = paste0(getter, collapse = "\n"))
# }

proc_output = function(res)
{
  if (length(res[[1]]) == 0) return(NULL)
  x = res[[2]]; out = list(res[[1]])
  if (any(diff(sapply(x, length)) != 0)) return(res)
  out = append(out, rw2cw(x))
  xnames = names(x[[1]])
  if (is.null(xnames) || length(xnames) != length(x[[1]]))
    xnames = paste0("X", 1:length(x[[1]]))
  names(out) = c("Time", xnames)
  attr(out, "row.names") = c(NA, -length(out[[1]]))
  class(out) = "data.frame"
  return(out)
}

rw2cw = function(x)
{
  lapply(lapply(1:length(x[[1]]),
                function(i) lapply(x, function(a) a[[i]])),
         unlist)
}

######################################################################################
######################################################################################
######################################################################################
######################################################################################
############## needs editing

##add parameter names that are vectors
# MACRO_VECTOR_PARS<- c("prec","epsilon","mu_A","theta","gamma_V",
#                       "phi_HV","phi_VH","mu_V_temp","mu_V_hum",
#                       "c_temp") 

MACRO_VECTOR_PARS<- c("prec","epsilon","mu_A","theta","gamma_V",
                      "phi_HV","phi_VH","mu_V_temp","mu_V_hum",
                      "c_temp","c_rain","mu_A_rain","a_rain","gam_V_temp","a_hum")



##add parameter names that are NOT vectors
MACRO_NONVECTOR_PARS<- c("ts","alpha","eta","gamH","sigH","a","f","K","delta")


JLmake_pars_decl = function(pars)
{

  ##build model text for C code

    res = paste0("double ", paste(MACRO_NONVECTOR_PARS, collapse = ", "), ";\n")

    for(par in MACRO_VECTOR_PARS){
        res = paste0(res, "vector<double> ",par,";\n")
    }

    for(par in MACRO_VECTOR_PARS){
        res = paste0(res, "double hw_",par,"[] = {",paste(as.numeric(unlist(model_input[par])),sep='',collapse=','), "};\n")
    }

    #can use this instead if getting C++ errors, easier to see error output on syntax in the console
    # for(rr in negNames)
    #   res = paste0(res, "double hw_indexp_m",rr,"[]= {1};\n")
    # res = paste0(res, "double hw_indexp_0[]= {1};\n")
    # for(rr in posNames)
    #   res = paste0(res, "double hw_indexp_p",rr,"[]= {1};\n")

    fileConn<-file("c_model_make_pars_decl.code.cpp"); writeLines(res, fileConn); close(fileConn)

    return(res)
}

JLmake_pars_setter = function(pars)
{

  ##build model text for C code

    ipars<- paste0("double ", MACRO_NONVECTOR_PARS, collapse = ", ")
    body<- paste0("odeintr::", MACRO_NONVECTOR_PARS, " = ", MACRO_NONVECTOR_PARS, collapse = ";\n")

    body<- paste0(body,";\n")

    for(par in MACRO_VECTOR_PARS){
        body<- paste0(body, "odeintr::",par,"=std::vector<double>(odeintr::hw_",par,", odeintr::hw_",par," + sizeof(odeintr::hw_",par,") / sizeof(double) );\n")
    }   

    code<- read_template("pars_setter_template")
    code<- sub("__PARS__", ipars, code)
    code<- sub("__BODY__", body, code)
    res<- paste0(code, collapse = "\n")

    fileConn<- file("c_model_make_pars_setter.code.cpp"); writeLines(res, fileConn); close(fileConn)

    return(res)
}

##JL
JLmake_pars_getter = function(pars)
{

  ##build model text for C code

    body<- paste0("out[\"", MACRO_NONVECTOR_PARS, "\"] = odeintr::", MACRO_NONVECTOR_PARS, collapse = ";\n")

    for(par in MACRO_VECTOR_PARS){
        body<- paste0(body, ";\n out[\"",par,"\"] = odeintr::",par,";");
    }   
 
    code<- read_template("pars_getter_template")
    code<- sub("__BODY__", body, code)
    res<- paste0(code, collapse = "\n")

    fileConn<-file("c_model_make_pars_getter.code.cpp"); writeLines(res, fileConn); close(fileConn)

    return(res)
}

### define the model 

print(paste0("Defining extra code C++...."))

SYSDIM<-13 ##number of equations

# C++ code for ODE system
SEIRsys = '

  \t /* int tt= int(t*1/ts); t comes in days, but step is ts, and parameters in scale of ts */

  /* t comes in days, but step is ts, and parameters in scale of ts */
  /* delta is the shift in time to consider the climate at a different time period */
  \t int tt= int((t+delta)*1/ts);  /*  note that delta comes negative */

  \t double SH= x[0];
  \t double EH= x[1];
  \t double IH= x[2];
  \t double RH= x[3];
  \t double AV= x[4];
  \t double SV= x[5];
  \t double EV= x[6];
  \t double IV= x[7]; 

  \t double INC= x[8];
  \t double INCCS= x[9];

  \t double R0= x[10];
  \t double VLS= x[11];
  \t double EIC= x[12];

  \t double NH= SH+EH+IH+RH;
  \t double NV= SV+EV+IV;

  \t double R= prec[tt];

  \t double phiVH= phi_VH[tt];
  \t double phiHV= phi_HV[tt];
  \t double cV= c_temp[tt]; 
  \t double thetaV= theta[tt];
  \t double epsA= epsilon[tt];
  \t double muA= eta * mu_A[tt] * pow(1 + mu_A_rain[tt],eta);
  \t double muV= eta * mu_V_temp[tt] * pow(1 + mu_V_hum[tt],eta); 
  \t double gamV= alpha * gamma_V[tt]; /* * pow(1 + gam_V_temp[tt], alpha); */
  \t double aV= a; /* * pow(1 + a_hum[tt],eta); */ 

  double NEAR_ZERO= 0.00000000001;
  double NEAR_ONE= 0.999999999;
  double MIN_MUV= 1/60; /* assume max life exp is 60 days */
  double MIN_MUA= 1/120; /* assume max life exp is 120 days */
  if(phiVH<0.0) phiVH=NEAR_ZERO;
  if(phiVH>1.0) phiVH=NEAR_ONE;
  if(phiHV<0.0) phiHV=NEAR_ZERO;
  if(phiHV>1.0) phiHV=NEAR_ONE;
  if(cV<0.0) cV=NEAR_ZERO;
  if(cV>1.0) cV=NEAR_ONE;
  if(epsA<0.0) epsA=NEAR_ZERO;
  if(muA<MIN_MUA) muA=MIN_MUA;
  \\\ if(muA>1.0) muA=NEAR_ONE;
  if(muV<MIN_MUV) muV=MIN_MUV;
  \\\ if(muV>1.0) muV=NEAR_ONE;
  if(gamV<0.0) muV=NEAR_ZERO;
  if(aV<0.0) aV=NEAR_ZERO;

  \t double lamVH = aV * phiVH * IV * SH / NH;
  \t double lamHV = aV * phiHV * IH * SV / NH;
  \t double V = SV + EV + IV;

  \t double r0 = (V/NH)*aV*aV*phiHV*phiVH*gamV*gamH/(muV*(sigH+0.0)*(gamH+0.0)*(gamV+muV));

  \t  dxdt[0] = -lamVH;            //SH
  \t  dxdt[1] = lamVH - gamH*EH;   //EH
  \t  dxdt[2] = gamH*EH -sigH*IH;   //IH
  \t  dxdt[3] = sigH*IH;            //RH

  \t  dxdt[4] = cV*f*thetaV*(1.0-AV/(K*(R+1)))*NV -epsA*AV -muA*AV ; //AV

  \t  dxdt[5] = epsA*AV -lamHV -muV*SV; ;    //SV
  \t  dxdt[6] = lamHV -gamV*EV -muV*EV;    //EV
  \t  dxdt[7] = gamV*EV -muV*IV;            //IV

  \t  dxdt[8] = gamH*EH - INC; //Inc
  \t  dxdt[9] = gamH*EH; //IncCS

  \t /* these next ones are here to save key variables that can be used later (e.g for likelihood) */
  \t  dxdt[10] = r0 - R0; //R0
  \t  dxdt[11] = 1.0/muV - VLS; //Vls
  \t  dxdt[12] = 1.0/gamV - EIC; //EIC

'

##set some extra code like headers and global functions
headers<- "#include <vector>\n
           #include <stdio.h>\n
           #include <stdlib.h>\n
           #include <fstream>\n
           #include <sstream>\n
           #include <iostream>\n
           #include <iterator>\n
           using namespace std;\n
           \n
           std::vector<double> &split(const string &s, char delim, std::vector<double> &elems) {\n
            \t stringstream ss(s);\n
            \t string item;\n
            \t while (std::getline(ss, item, delim)) elems.push_back(atof(item.c_str()));\n
            \t return elems;
            }\n
            "

##compile and create the actual functions to call the model
##need doing only once per MCMC (not each step!)
print(paste0("Now compiling C++...."))

thesepars<- c(MACRO_VECTOR_PARS,MACRO_NONVECTOR_PARS)

# code<- mycompile_sys(name="ODEmodel", sys=SEIRsys,
#                      pars=thesepars,
#                      headers=headers,
#                      method='rk4', sys_dim=SYSDIM, compile=TRUE)

code<- mycompile_sys(name="ODEmodel", sys=SEIRsys,
                     pars=thesepars,
                     headers=headers,
                     method='rk4', sys_dim=SYSDIM)

print(paste0("Code saved...."))
print(paste0("Code saved...."))

fileConn<-file("c_model.code.cpp")
writeLines(code, fileConn)
close(fileConn)
print(paste0("Done!"))
