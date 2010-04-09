libobjs = switch_option.o pswitch.o param.o autohelp.o autoversion.o

switch_option.o : switch_option.hh
pswitch.o : pswitch.hh switch_option.hh
param.o : param.hh pswitch.hh switch_option.hh
autohelp.o : autohelp.hh switch_option.hh pswitch.hh
autoversion.o : autoversion.hh switch_option.hh pswitch.hh

libparammm.a: $(libobjs)
	ar ruv libparammm.a $?
	ranlib libparammm.a

test : test.o libparammm.a
	$(CXX) -o test test.o -L. -lparammm
