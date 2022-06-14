/*
  Copyright (C) 2017 Julian Späth
  ----------------------------------------------------------------------------
  This file is part of AMMSoCK.

  AMMSoCK is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version. AMMSoCK is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
  details. You should have received a copy of the GNU General Public License
  along with AMMSoCK. If not, see <https://www.gnu.org/licenses/>.
*/

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "Python.h"
#include "ammsockNLP.hpp"
#include "numpy/arrayobject.h"

using namespace Ipopt;
static PyObject *Reduction(PyObject *self, PyObject *args) {
  // initialize
  Py_Initialize();

  SmartPtr<ammsockNLP> mynlp = new ammsockNLP();

  //----------------------------------------------------------------------------------------------------------//
  // check number of arguments
  if (PyTuple_Size(args) != 3) {
    std::cerr << "performReduction - Error: Function call is "
                 "performReduction(rpv,init,param)."
              << std::endl;
  }

  PyObject *rpv_ob;
  PyObject *init_ob;
  PyObject *param_ob;

  //----------------------------------------------------------------------------------------------------------//
  // Check for arguments
  if (!PyArg_ParseTuple(args, "OOO", &rpv_ob, &init_ob, &param_ob)) {
    std::cerr << "performReduction - Error: wrong format." << std::endl;
    exit(0);
  }

  // Check for list in first argument
  if (!PyList_Check(rpv_ob)) {
    std::cerr << "performReduction - Error: first argument need to be a list."
              << std::endl;
    exit(0);
  }

  // check for list in second argument
  if (!PyList_Check(init_ob)) {
    std::cerr << "performReduction - Error: second argument need to be a list."
              << std::endl;
    exit(0);
  }

  // Check for dictionary in third argument
  if (!PyDict_Check(param_ob)) {
    std::cerr
        << "performReduction - Error: third argument need to be a dictionary."
        << std::endl;
    exit(0);
  }

  //----------------------------------------------------------------------------------------------------------//
  // Check if arguments are empty

  if (PyList_Size(rpv_ob) == 0) {
    std::cerr << "performReduction - Error: rpv has to be non-empty."
              << std::endl;
    exit(0);
  }

  if (PyList_Size(init_ob) == 0) {
    std::cerr << "performReduction - Error: init has to be non-empty."
              << std::endl;
    exit(0);
  }

  if (PyDict_Size(param_ob) == 0) {
    std::cerr << "performReduction - Error: param has to be non-empty."
              << std::endl;
    exit(0);
  }

  //----------------------------------------------------------------------------------------------------------//
  // transfer from python dictionary "param" to c++

  Py_ssize_t pos;
  PyObject *key;
  char *str;

  // hcoll//
  // check if hcoll is empty
  pos = 0;
  str = "hcoll";
  key = Py_BuildValue("s", str);
  if (PyDict_Contains(param_ob, key) == 0) {
    std::cerr << "performReduction - Error: hcoll field is empty." << std::endl;
    exit(0);
  }

  // transfer
  PyObject *hcoll_ = PyDict_GetItemString(param_ob, str);
  double hcoll;

  if (PyLong_Check(hcoll_) || PyFloat_Check(hcoll_)) {
    hcoll = PyFloat_AsDouble(hcoll_);
  } else {
    std::cerr << "Error in Python dictionary: hcoll has to be a double"
              << std::endl;
    exit(0);
  }

  // hfixed//
  // check if scale is empty
  pos = 0;
  str = "hfixed";
  key = Py_BuildValue("s", str);
  if (PyDict_Contains(param_ob, key) == 0) {
    std::cerr << "performReduction - Error: hixed field is empty." << std::endl;
    exit(0);
  }

  // transfer

  PyObject *hfixed_ = PyDict_GetItemString(param_ob, str);
  double hfixed;

  if (PyLong_Check(hfixed_) || PyFloat_Check(hfixed_)) {
    hfixed = PyFloat_AsDouble(hfixed_);
  } else {
    std::cerr << "Error in Python dictionary: hfixed has to be a double"
              << std::endl;
    exit(0);
  }

  // scale//
  // check if scale is empty
  pos = 0;
  str = "scale";
  key = Py_BuildValue("s", str);
  if (PyDict_Contains(param_ob, key) == 0) {
    std::cerr << "performReduction - Error: scale field is empty." << std::endl;
    exit(0);
  }

  // transfer

  PyObject *scale_ = PyDict_GetItemString(param_ob, str);
  double scale;

  if (PyLong_Check(scale_) || PyFloat_Check(scale_)) {
    scale = PyFloat_AsDouble(scale_);
  } else {
    std::cerr << "Error in Python dictionary: scale has to be a double"
              << std::endl;
    exit(0);
  }

  // atomCons//
  // check if atomCons is empty
  pos = 0;
  str = "atomCons";
  key = Py_BuildValue("s", str);
  if (PyDict_Contains(param_ob, key) == 0) {
    std::cerr << "performReduction - Error: atomCons field is empty."
              << std::endl;
    exit(0);
  }

  // transfer to c++ array
  PyObject *list = PyList_New(0);
  PyObject *atomCons_ = PyDict_GetItemString(param_ob, str);

  if (PyLong_Check(atomCons_) || PyFloat_Check(atomCons_)) {
    PyList_Append(list, atomCons_);
  } else if (PyList_Check(atomCons_)) {
    list = atomCons_;
  } else {
    std::cerr
        << "Error in Python dictionary: atomCons has to be a list or a double"
        << std::endl;
    exit(0);
  }

  double atomCons[int(PyList_Size(list))];

  for (int k = 0; k < int(PyList_Size(list)); k++) {
    atomCons[k] = PyFloat_AsDouble(PyList_GetItem(list, pos));
    pos++;
  }

  // maxit//
  // check if maxit is empty and transfer to c++ double
  int maxit;
  str = "maxit";
  key = Py_BuildValue("s", str);
  if (PyDict_Contains(param_ob, key) == 0) {
    std::cerr << "performReduction - Error: maxit field is empty." << std::endl;
    ;
    exit(0);
  } else {
    maxit = PyLong_AsDouble(PyDict_GetItemString(param_ob, str));
  }

  // exactHessian//
  // check if ecaxt is empty and transfer to c++ double
  int exactHessian;
  str = "exactHessian";
  key = Py_BuildValue("s", str);
  if (PyDict_Contains(param_ob, key) == 0) {
    std::cerr << "performReduction - Error: exactHessian field is empty."
              << std::endl;
  } else {
    exactHessian = PyLong_AsDouble(PyDict_GetItemString(param_ob, str));
  }

  //----------------------------------------------------------------------------------------------------------//
  // transfer from python list rpv to c++ array

  double rpv_tmp[int(PyList_Size(rpv_ob))];
  double *rpv;
  pos = 0;

  for (int k = 0; k < int(PyList_Size(rpv_ob)); k++) {
    rpv_tmp[k] = PyFloat_AsDouble(PyList_GetItem(rpv_ob, pos));
    pos++;
  }
  rpv = rpv_tmp;

  //----------------------------------------------------------------------------------------------------------//
  // transfer from python list init to c++ array
  double init_tmp[int(PyList_Size(init_ob))];
  double *init;
  pos = 0;

  for (int k = 0; k < int(PyList_Size(init_ob)); k++) {
    init_tmp[k] = PyFloat_AsDouble(PyList_GetItem(init_ob, pos));
    pos++;
  }
  init = init_tmp;

  //----------------------------------------------------------------------------------------------------------//

  mynlp->setParameter(rpv, init, hcoll, hfixed, scale, atomCons);
  SmartPtr<IpoptApplication> app = new IpoptApplication();

  ApplicationReturnStatus status;
  status = app->Initialize();
  app->Options()->SetStringValue("mu_strategy", "monotone");

  if (exactHessian == 0) {
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
  } else {
    app->Options()->SetStringValue("hessian_approximation", "exact");
  }

  app->Options()->SetStringValue("warm_start_init_point", "yes");
  app->Options()->SetStringValue("expect_infeasible_problem", "yes");
  app->Options()->SetIntegerValue("max_iter", 10000);
  app->Options()->SetIntegerValue("print_level", 0);

  if (status != Solve_Succeeded) {
    std::cerr << std::endl
              << std::endl
              << "*** Error during initialization!" << std::endl;
  }

  int it = 1;
  status = app->OptimizeTNLP(mynlp);
  while ((status != Solve_Succeeded) &&
         (status != Solved_To_Acceptable_Level) && (it <= maxit)) {
    status = app->ReOptimizeTNLP(mynlp);
    it++;
  }

  double result[NOP + 1];

  mynlp->getSolution(result);

  //---return Array

  // npy_intp dims[1] = {NOP+1};
  // PyObject *ret = PyArray_SimpleNew(1,dims,NPY_DOUBLE);
  // memcpy(PyArray_DATA(ret),result,sizeof(result));
  //---------------

  //---return List

  PyObject *ret = PyList_New(NOP + 1);
  pos = 0;
  for (int k = 0; k < NOP + 1; k++) {
    PyList_SetItem(ret, pos, PyFloat_FromDouble(result[pos]));
    pos++;
  }
  //---------------

  return Py_BuildValue("iiO", status, it, ret);
}

static PyMethodDef SpamMethods[] = {
    {"Reduction", Reduction, METH_VARARGS, "Execute a shell command."},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef performReduction = {
    PyModuleDef_HEAD_INIT, "performReduction", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module,
             or -1 if the module keeps state in global variables. */
    SpamMethods};

PyMODINIT_FUNC PyInit_performReduction(void) {
  return PyModule_Create(&performReduction);
}
