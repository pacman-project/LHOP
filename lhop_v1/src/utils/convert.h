/********************************************************************
	File:		convert.h
	Created:	2004/02/24
	Author:		Daniel Eloff
	
	Purpose:	Functions for converting string to integer/floating-point types

	License:	This code is my work. Use it freely
				so long as you give me credit for it.
*********************************************************************/


#ifndef CONVERT_INCLUDED
#define CONVERT_INCLUDED

#include <limits>
#include <math.h>

class InvalidConversionException {};

class Convert
{
private:
	typedef unsigned long size_type;

	static double __ToDouble(const char * str, size_type size, size_type& i)
	{
		unsigned long num = 0;
		unsigned long decimals = 0;
		double decimalPlaces;

		if(size == 0) //don't read beyond the end of the number!
			return 0;

		//Check if this is a negative number
		int sign = 1;
		if(str[i] == '-')
		{
			sign = -1;
			++i; //scroll past sign
		}
		else if(str[i] == '+')
			++i; //scrol pas the unescessary sign

		//Convert each digit char and add into result.
		for ( ; i < size && str[i] >= '0' && str[i] <= '9'; ++i) 
			num = (num * 10) + (str[i] - '0');

		//Do the same for the decimals
		if(i < size && str[i] == '.')
		{
			decimalPlaces = i; //just storing this temporarily
			for(++i; i < size && str[i] >= '0' && str[i] <= '9'; ++i) 
				decimals = (decimals * 10) + (str[i] - '0');
            decimalPlaces = pow(10.0,int(i-decimalPlaces-1)); //get the number of decimal places
			return (num+decimals/decimalPlaces)*sign;
		}
		return num*static_cast<double>(sign);
	}

public:
	//Convert string to a double
	static double ToDouble(const char * str, size_type size)
	{
		size_type i = 0;
		double exp = 0;

		//scroll past whitespace/characters until we find the start of a number
		for(; i < size; str++, ++i)
		{
			//these conditions are ordered for effeciency
			if(*str <= '9' && *str >= '+') //first condition excludes all ASCII chars above '9', second all chars under '+'
			{
				if(*str >= '0') //definately got a number
					break;
				else if((*str == '-' || *str == '+') && i+1 < size)
				{
					if(str[1] <= '9' && str[1] >= '0')
						break; //got a signed number
					else if(str[1] == '.' && i+2 < size)
					{
						if(str[2] <= '9' && str[2] >= '0')
							break; //got a decimal number
					}
				}
				else if(*str == '.')
				{
					if(str[1] <= '9' && str[1] >= '0')
						break; //got a decimal number
				}
			}
		}
		size -= i; //adjust size
		i=0; //reset i to zero

		//Get the base
		double base = __ToDouble(str,size,i);

		//Look for exponents
		if(i < size && str[i] == 'e' || str[i] == 'E')
			exp = __ToDouble(str,size,++i);

		if(exp != 0)
			return base*pow(double(10),exp);
		else if(i == 0)
			throw InvalidConversionException();
		return base;
	}

	//Convet a string to an integral type
	template<class I>
	static I ToInteger(const char * str, unsigned long size)
	{
		size_type i = 0;
		if(size == 0)
			throw InvalidConversionException();
		I num = 0;

		//scroll past whitespace/characters until we find the start of a number
		for(; i < size; str++, ++i)
		{
			//these conditions are ordered for effeciency
			if(*str <= '9' && *str >= '+') //first condition excludes all ASCII chars above '9', second all chars under '+'
			{
				if(*str >= '0') //definately got a number
					break;
				else if((*str == '-' || *str == '+') && i+1 < size)
				{
					if(str[1] <= '9' && str[1] >= '0')
						break; //got a signed number
				}
			}
		}
		size -= i; //adjust size
		i=0; //reset i to zero

		//Check if this is a negative number
		int sign = 1;
		if(*str == '-')
		{
			if(std::numeric_limits<I>::is_signed)
				sign = -1; //don't make it a negative number unless we're converting to a signed type
			++i; //scroll past sign
		}
		else if(*str == '+')
			++i; //scroll past unescessary sign

		//Convert each digit char and add into result.
		for ( ; i < size && str[i] >= '0' && str[i] <='9'; ++i) 
			num = (num * 10) + (str[i] - '0');

		if(i == 0)
			throw InvalidConversionException();
		return num*sign;
	}
};

#endif
