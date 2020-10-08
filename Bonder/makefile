boostH=/usr/include/boost
boostLib=/usr/lib/x86_64-linux-gnu

bonder : Analize.cpp fill.cpp iface.cpp LLi.cpp main.cpp readwfn.cpp stdafx.cpp output.cpp
	g++ -o bond Analize.cpp fill.cpp iface.cpp LLi.cpp main.cpp readwfn.cpp stdafx.cpp output.cpp -I . -I $(boostH)   -std=c++11  $(boostLib)/libboost_thread.a $(boostLib)/libboost_system.a $(boostLib)/libboost_chrono.a -pthread  -Ofast #-fsanitize=undefined
