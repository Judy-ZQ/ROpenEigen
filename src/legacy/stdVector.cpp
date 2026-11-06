#include "ROpenEigen_types.h"


RCPP_MODULE(stdVectorModule){
  class_<myStdVector<int>>("stdVectorOfi")
  .default_constructor()
  .factory<SEXP>(newstdVector<int>)
  // .constructor<int>()
  // .constructor<myStdVector<int>>()
  // exposing member functions
  .method("size", &myStdVector<int>::size)
  .method("max_size", &myStdVector<int>::max_size)
  .method("capacity", &myStdVector<int>::capacity)
  .method("empty", &myStdVector<int>::empty)
  .method("reserve", &myStdVector<int>::reserve)
  .method("pop_back", &myStdVector<int>::pop_back)
  .method("clear", &myStdVector<int>::clear)
  // .method("XPtr_",&myStdVector<int>::getXPtr)
     .method("XPtr_",&stdVector_getXPtr<int>)
     .method("outToConsole",&stdVectorOutToConsole<int>)
  
  // exposing const member functions
  .const_method("back", &myStdVector<int>::back)
  .const_method("front", &myStdVector<int>::front)
  .const_method("at", &myStdVector<int>::at1 )
  
  // exposing free functions taking a
  // std::vector<double>* as their first
  // argument
  .method("assign", &stdVector_assign<int>)
  .method("insert", &stdVector_insert<int>)
  .method("push_back", &stdVector_push_back<int>)
  .method("resize", &stdVector_resize<int>)
  .method("std2r", &stdVector_std2r<int>)
  
  // special methods for indexing
  .const_method("[[", &myStdVector<int>::at1)
  .method("[[<-", &stdVector_set<int>)
  
  ;
  
  class_<myStdVector<float>>("stdVectorOff")
    .default_constructor()
    .factory<SEXP>(newstdVector<float>)
  // .constructor<float>()
  // .constructor<myStdVector<float>>()
  // exposing member functions
  .method("size", &myStdVector<float>::size)
  .method("max_size", &myStdVector<float>::max_size)
  .method("capacity", &myStdVector<float>::capacity)
  .method("empty", &myStdVector<float>::empty)
  .method("reserve", &myStdVector<float>::reserve)
  .method("pop_back", &myStdVector<float>::pop_back)
  .method("clear", &myStdVector<float>::clear)
  // .method("XPtr_",&myStdVector<float>::getXPtr)
     .method("XPtr_",&stdVector_getXPtr<float>)
     .method("outToConsole",&stdVectorOutToConsole<float>)
  
  // exposing const member functions
  .const_method("back", &myStdVector<float>::back)
  .const_method("front", &myStdVector<float>::front)
  .const_method("at", &myStdVector<float>::at1 )
  
  // exposing free functions taking a
  // std::vector<float>* as their first
  // argument
  .method("assign", &stdVector_assign<float>)
  .method("insert", &stdVector_insert<float>)
  .method("push_back", &stdVector_push_back<float>)
  .method("resize", &stdVector_resize<float>)
  .method("std2r", &stdVector_std2r<float>)
  
  // special methods for indexing
  .const_method("[[", &myStdVector<float>::at1)
  .method("[[<-", &stdVector_set<float>)
  
  ;
  
  class_<myStdVector<double>>("stdVectorOfd")
    .default_constructor()
    .factory<SEXP>(newstdVector<double>)
  // .constructor<double>()
  // .constructor<myStdVector<double>>()
  // exposing member functions
  .method("size", &myStdVector<double>::size)
  .method("max_size", &myStdVector<double>::max_size)
  .method("capacity", &myStdVector<double>::capacity)
  .method("empty", &myStdVector<double>::empty)
  .method("reserve", &myStdVector<double>::reserve)
  .method("pop_back", &myStdVector<double>::pop_back)
  .method("clear", &myStdVector<double>::clear)
  // .method("XPtr_",&myStdVector<double>::getXPtr)
     .method("XPtr_",&stdVector_getXPtr<double>)
     .method("outToConsole",&stdVectorOutToConsole<double>)
  
  // exposing const member functions
  .const_method("back", &myStdVector<double>::back)
  .const_method("front", &myStdVector<double>::front)
  .const_method("at", &myStdVector<double>::at1 )
  
  // exposing free functions taking a
  // std::vector<double>* as their first
  // argument
  .method("assign", &stdVector_assign<double>)
  .method("insert", &stdVector_insert<double>)
  .method("push_back", &stdVector_push_back<double>)
  .method("resize", &stdVector_resize<double>)
  .method("std2r", &stdVector_std2r<double>)
  
  // special methods for indexing
  .const_method("[[", &myStdVector<double>::at1)
  .method("[[<-", &stdVector_set<double>)
  
  ;
  
  class_<myStdVector<char>>("stdVectorOfc")
    .default_constructor()
    .factory<SEXP>(newstdVector<char>)
  // .constructor<char>()
  // .constructor<myStdVector<char>>()
  // exposing member functions
  .method("size", &myStdVector<char>::size)
  .method("max_size", &myStdVector<char>::max_size)
  .method("capacity", &myStdVector<char>::capacity)
  .method("empty", &myStdVector<char>::empty)
  .method("reserve", &myStdVector<char>::reserve)
  .method("pop_back", &myStdVector<char>::pop_back)
  .method("clear", &myStdVector<char>::clear)
  // .method("XPtr_",&myStdVector<char>::getXPtr)
     .method("XPtr_",&stdVector_getXPtr<char>)
     .method("outToConsole",&stdVectorOutToConsole<char>)
  
  // exposing const member functions
  .const_method("back", &myStdVector<char>::back)
  .const_method("front", &myStdVector<char>::front)
  .const_method("at", &myStdVector<char>::at1 )
  
  // exposing free functions taking a
  // std::vector<double>* as their first
  // argument
  .method("assign", &stdVector_assign<char>)
  .method("insert", &stdVector_insert<char>)
  .method("push_back", &stdVector_push_back<char>)
  .method("resize", &stdVector_resize<char>)
  .method("std2r", &stdVector_std2r<char>)
  
  // special methods for indexing
  .const_method("[[", &myStdVector<char>::at1)
  .method("[[<-", &stdVector_set<char>)
  
  ;
  
  class_<myStdVector<unsigned char>>("stdVectorOfuc")
    .default_constructor()
    .factory<SEXP>(newstdVector<unsigned char>)
  // .constructor<unsigned char>()
  // .constructor<myStdVector<unsigned char>>()
  // exposing member functions
  .method("size", &myStdVector<unsigned char>::size)
  .method("max_size", &myStdVector<unsigned char>::max_size)
  .method("capacity", &myStdVector<unsigned char>::capacity)
  .method("empty", &myStdVector<unsigned char>::empty)
  .method("reserve", &myStdVector<unsigned char>::reserve)
  .method("pop_back", &myStdVector<unsigned char>::pop_back)
  .method("clear", &myStdVector<unsigned char>::clear)
  // .method("XPtr_",&myStdVector<unsigned char>::getXPtr)
     .method("XPtr_",&stdVector_getXPtr<unsigned char>)
     .method("outToConsole",&stdVectorOutToConsole<unsigned char>)
  
  // exposing const member functions
  .const_method("back", &myStdVector<unsigned char>::back)
  .const_method("front", &myStdVector<unsigned char>::front)
  .const_method("at", &myStdVector<unsigned char>::at1 )
  
  // exposing free functions taking a
  // std::vector<double>* as their first
  // argument
  .method("assign", &stdVector_assign<unsigned char>)
  .method("insert", &stdVector_insert<unsigned char>)
  .method("push_back", &stdVector_push_back<unsigned char>)
  .method("resize", &stdVector_resize<unsigned char>)
  .method("std2r", &stdVector_std2r<unsigned char>)
  
  // special methods for indexing
  .const_method("[[", &myStdVector<unsigned char>::at1)
  .method("[[<-", &stdVector_set<unsigned char>)
  
  ;
  
  class_<myStdVector<myStdVector<char>>>("stdVectorOfVectorOfc")
    .default_constructor()
    .factory<SEXP>(newstdVector<myStdVector<char>>)
    .method("size", &myStdVector<myStdVector<char>>::size)
    .method("max_size", &myStdVector<myStdVector<char>>::max_size)
    .method("capacity", &myStdVector<myStdVector<char>>::capacity)
    .method("empty", &myStdVector<myStdVector<char>>::empty)
    .method("reserve", &myStdVector<myStdVector<char>>::reserve)
    .method("pop_back", &myStdVector<myStdVector<char>>::pop_back)
    .method("clear", &myStdVector<myStdVector<char>>::clear)
  // .method("XPtr_",&myStdVector<myStdVector<char>>::getXPtr)
     .method("XPtr_",&stdVector_getXPtr<myStdVector<char>>)
     .const_method("back", &myStdVector<myStdVector<char>>::back)
     .const_method("front", &myStdVector<myStdVector<char>>::front)
     .const_method("at", &myStdVector<myStdVector<char>>::at1 )
     .method("assign", &stdVector_assign<myStdVector<char>>)
     .method("insert", &stdVector_insert<myStdVector<char>>)
     .method("push_back", &stdVector_push_back<myStdVector<char>>)
     .method("resize", &stdVector_resize<myStdVector<char>>)
     .const_method("[[", &myStdVector<myStdVector<char>>::at1)
     .method("[[<-", &stdVector_set<myStdVector<char>>)
  ;
  
  class_<myStdVector<std::string>>("stdVectorOfString")
    .default_constructor()
    .factory<SEXP>(newstdVector<std::string>)
  // exposing member functions
  .method("size", &myStdVector<std::string>::size)
  .method("max_size", &myStdVector<std::string>::max_size)
  .method("capacity", &myStdVector<std::string>::capacity)
  .method("empty", &myStdVector<std::string>::empty)
  .method("reserve", &myStdVector<std::string>::reserve)
  .method("pop_back", &myStdVector<std::string>::pop_back)
  .method("clear", &myStdVector<std::string>::clear)
  // .method("XPtr_",&myStdVector<std::string>::getXPtr)
     .method("XPtr_",&stdVector_getXPtr<std::string>)
     .method("outToConsole",&stdVectorOutToConsole<std::string>)
  
  // exposing const member functions
  .const_method("back", &myStdVector<std::string>::back)
  .const_method("front", &myStdVector<std::string>::front)
  .const_method("at", &myStdVector<std::string>::at1 )
  
  // exposing free functions taking a
  // std::vector<double>* as their first
  // argument
  // .method("assign", &stdVector_assign<std::string>) //不能添加？提示=指代不明确
  // .method("insert", &stdVector_insert<cv::String>) //不能添加？
     .method("push_back", &stdVector_push_back<std::string>)
     .method("resize", &stdVector_resize<std::string>)
     .method("std2r", &stdVector_std2r<std::string>)
  
  // special methods for indexing
  .const_method("[[", &myStdVector<std::string>::at1)
  .method("[[<-", &stdVector_set<std::string>)
  
  ;
}
