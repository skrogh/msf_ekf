//
// Minimal library with some concurrent queues and buffers.
// Compile and link with Boost (circular_buffer is headder only. Thread is not)
//

#pragma once
#include <queue>
#include <boost/thread.hpp>
#include <boost/circular_buffer.hpp>


template<typename Data>
class ConcurrentQueue
{
private:
	std::queue<Data> the_queue;
	mutable boost::mutex the_mutex;
	boost::condition_variable the_condition_variable;
public:
	ConcurrentQueue(){}

	ConcurrentQueue<Data>& operator=(const ConcurrentQueue<Data>& other)
	{
		boost::mutex::scoped_lock lock_other(other.the_mutex);
		boost::mutex::scoped_lock lock(the_mutex);
		the_queue = other.the_queue;
		lock.unlock();
		lock_other.unlock();
		the_condition_variable.notify_one();
	}
	
	ConcurrentQueue(const ConcurrentQueue<Data>& other)
	{
		boost::mutex::scoped_lock lock_other(other.the_mutex);
		boost::mutex::scoped_lock lock(the_mutex);
		the_queue = other.the_queue;
		lock.unlock();
		lock_other.unlock();
		the_condition_variable.notify_one();
	}

	void push(Data const& data)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		the_queue.push(data);
		lock.unlock();
		the_condition_variable.notify_one();
	}

	bool empty() const
	{
		boost::mutex::scoped_lock lock(the_mutex);
		return the_queue.empty();
	}

	bool peek(Data& peeked_value) const
	{
		boost::mutex::scoped_lock lock(the_mutex);
		if(the_queue.empty())
		{
			return false;
		}
		
		peeked_value=the_queue.front();
		return true;
	}

	bool popNonBlocking(Data& popped_value)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		if(the_queue.empty())
		{
			return false;
		}
		
		popped_value=the_queue.front();
		the_queue.pop();
		return true;
	}

	void blockWhileEmpty()
	{
		boost::mutex::scoped_lock lock(the_mutex);
		while(the_queue.empty())
		{
			the_condition_variable.wait(lock);
		}
	}

	void popBlocking(Data& popped_value)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		while(the_queue.empty())
		{
			the_condition_variable.wait(lock);
		}
		
		popped_value=the_queue.front();
		the_queue.pop();
	}

};


template<typename Data,int Length>
class ConcurrentFifo
{
private:
	boost::circular_buffer<Data> the_queue;
	mutable boost::mutex the_mutex;
	boost::condition_variable the_condition_variable;
public:
	ConcurrentFifo(): the_queue(Length)
	{
		
	}
	
	ConcurrentFifo<Data,Length>& operator=(const ConcurrentFifo<Data,Length>& other)
	{
		boost::mutex::scoped_lock lock_other(other.the_mutex);
		boost::mutex::scoped_lock lock(the_mutex);
		the_queue = other.the_queue;
		lock.unlock();
		lock_other.unlock();
		the_condition_variable.notify_one();
	}
	
	ConcurrentFifo(const ConcurrentFifo<Data,Length>& other)
	{
		boost::mutex::scoped_lock lock_other(other.the_mutex);
		boost::mutex::scoped_lock lock(the_mutex);
		the_queue = other.the_queue;
		lock.unlock();
		lock_other.unlock();
		the_condition_variable.notify_one();
	}

	void push(Data const& data)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		the_queue.push_back(data);
		lock.unlock();
		the_condition_variable.notify_one();
	}

	bool empty() const
	{
		boost::mutex::scoped_lock lock(the_mutex);
		return the_queue.empty();
	}

	bool peek(Data& peeked_value) const
	{
		boost::mutex::scoped_lock lock(the_mutex);
		if(the_queue.empty())
		{
			return false;
		}
		
		peeked_value=the_queue.front();
		return true;
	}

	bool popNonBlocking(Data& popped_value)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		if(the_queue.empty())
		{
			return false;
		}
		
		popped_value=the_queue.front();
		the_queue.pop_front();
		return true;
	}

	void blockWhileEmpty()
	{
		boost::mutex::scoped_lock lock(the_mutex);
		while(the_queue.empty())
		{
			the_condition_variable.wait(lock);
		}
	}

	void popBlocking(Data& popped_value)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		while(the_queue.empty())
		{
			the_condition_variable.wait(lock);
		}
		
		popped_value=the_queue.front();
		the_queue.pop_front();
	}

};