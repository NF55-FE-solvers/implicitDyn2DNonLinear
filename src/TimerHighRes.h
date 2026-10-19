#ifndef TIMER_HIGH_RES_H
#define TIMER_HIGH_RES_H

#include <chrono>

class Timer
{
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start {};
    std::chrono::time_point<std::chrono::high_resolution_clock> m_stop {};

public:
    Timer() = default;

    void start(){ m_start = std::chrono::high_resolution_clock::now(); }
    void stop(){ m_stop = std::chrono::high_resolution_clock::now(); }

    auto getDurationMicro() const{
        std::chrono::duration<double,std::milli> d {(m_stop-m_start)*1000};
        return d.count();
    }    
};

#endif