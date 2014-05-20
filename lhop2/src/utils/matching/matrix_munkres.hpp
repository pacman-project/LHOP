/*
 *   Copyright (c) 2007 John Weaver
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 */

#include "matrix_munkres.h"

#include <cassert>
#include <cstdlib>
#include <algorithm>

/*export*/ template <class T>
MatrixMunkres<T>::MatrixMunkres() {
	m_rows = 0;
	m_columns = 0;
	m_matrix = nullptr;
}

/*export*/ template <class T>
MatrixMunkres<T>::MatrixMunkres(const MatrixMunkres<T> &other) {
	if ( other.m_matrix != nullptr ) {
		// copy arrays
		m_matrix = nullptr;
		resize(other.m_rows, other.m_columns);
		for ( int i = 0 ; i < m_rows ; i++ )
			for ( int j = 0 ; j < m_columns ; j++ )
				m_matrix[i][j] = other.m_matrix[i][j];
	} else {
		m_matrix = nullptr;
		m_rows = 0;
		m_columns = 0;
	}
}

/*export*/ template <class T>
MatrixMunkres<T>::MatrixMunkres(int rows, int columns) {
	m_matrix = nullptr;
	resize(rows, columns);
}

/*export*/ template <class T>
MatrixMunkres<T> &
MatrixMunkres<T>::operator= (const MatrixMunkres<T> &other) {
	if ( other.m_matrix != nullptr ) {
		// copy arrays
		resize(other.m_rows, other.m_columns);
		for ( int i = 0 ; i < m_rows ; i++ )
			for ( int j = 0 ; j < m_columns ; j++ )
				m_matrix[i][j] = other.m_matrix[i][j];
	} else {
		// free arrays
		for ( int i = 0 ; i < m_columns ; i++ )
			delete [] m_matrix[i];

		delete [] m_matrix;

		m_matrix = nullptr;
		m_rows = 0;
		m_columns = 0;
	}
	
	return *this;
}

/*export*/ template <class T>
MatrixMunkres<T>::~MatrixMunkres() {
	if ( m_matrix != nullptr ) {
		// free arrays
		for ( int i = 0 ; i < m_rows ; i++ )
			delete [] m_matrix[i];

		delete [] m_matrix;
	}
	m_matrix = nullptr;
}

/*export*/ template <class T>
void
MatrixMunkres<T>::resize(int rows, int columns) {
	if ( m_matrix == nullptr ) {
		// alloc arrays
		m_matrix = new T*[rows]; // rows
		for ( int i = 0 ; i < rows ; i++ )
			m_matrix[i] = new T[columns]; // columns

		m_rows = rows;
		m_columns = columns;
		clear();
	} else {
		// save array pointer
		T **new_matrix;
		// alloc new arrays
		new_matrix = new T*[rows]; // rows
		for ( int i = 0 ; i < rows ; i++ ) {
			new_matrix[i] = new T[columns]; // columns
			for ( int j = 0 ; j < columns ; j++ )
				new_matrix[i][j] = 0;
		}

		// copy data from saved pointer to new arrays
		int minrows = std::min<int>(rows, m_rows);
		int mincols = std::min<int>(columns, m_columns);
		for ( int x = 0 ; x < minrows ; x++ )
			for ( int y = 0 ; y < mincols ; y++ )
				new_matrix[x][y] = m_matrix[x][y];

		// delete old arrays
		if ( m_matrix != nullptr ) {
			for ( int i = 0 ; i < m_rows ; i++ )
				delete [] m_matrix[i];

			delete [] m_matrix;
		}

		m_matrix = new_matrix;
	}

	m_rows = rows;
	m_columns = columns;
}

/*export*/ template <class T>
void
MatrixMunkres<T>::identity() {
	assert( m_matrix != nullptr );

	clear();

	int x = std::min<int>(m_rows, m_columns);
	for ( int i = 0 ; i < x ; i++ )
		m_matrix[i][i] = 1;
}

/*export*/ template <class T>
void
MatrixMunkres<T>::clear() {
	assert( m_matrix != nullptr );

	for ( int i = 0 ; i < m_rows ; i++ )
		for ( int j = 0 ; j < m_columns ; j++ )
			m_matrix[i][j] = 0;
}

/*export*/ template <class T>
T 
MatrixMunkres<T>::trace() {
	assert( m_matrix != nullptr );

	T value = 0;

	int x = std::min<int>(m_rows, m_columns);
	for ( int i = 0 ; i < x ; i++ )
		value += m_matrix[i][i];

	return value;
}

/*export*/ template <class T>
MatrixMunkres<T>& 
MatrixMunkres<T>::transpose() {
	assert( m_rows > 0 );
	assert( m_columns > 0 );

	int new_rows = m_columns;
	int new_columns = m_rows;

	if ( m_rows != m_columns ) {
		// expand matrix
		int m = std::max<int>(m_rows, m_columns);
		resize(m,m);
	}

	for ( int i = 0 ; i < m_rows ; i++ ) {
		for ( int j = i+1 ; j < m_columns ; j++ ) {
			T tmp = m_matrix[i][j];
			m_matrix[i][j] = m_matrix[j][i];
			m_matrix[j][i] = tmp;
		}
	}

	if ( new_columns != new_rows ) {
		// trim off excess.
		resize(new_rows, new_columns);
	}

	return *this;
}

/*export*/ template <class T>
MatrixMunkres<T> 
MatrixMunkres<T>::product(MatrixMunkres<T> &other) {
	assert( m_matrix != nullptr );
	assert( other.m_matrix != nullptr );
	assert ( m_columns == other.m_rows );

	MatrixMunkres<T> out(m_rows, other.m_columns);

	for ( int i = 0 ; i < out.m_rows ; i++ ) {
		for ( int j = 0 ; j < out.m_columns ; j++ ) {
			for ( int x = 0 ; x < m_columns ; x++ ) {
				out(i,j) += m_matrix[i][x] * other.m_matrix[x][j];
			}
		}
	}

	return out;
}

/*export*/ template <class T>
T&
MatrixMunkres<T>::operator ()(int x, int y) {
	assert ( x >= 0 );
	assert ( y >= 0 );
	assert ( x < m_rows );
	assert ( y < m_columns );
	assert ( m_matrix != nullptr );
	return m_matrix[x][y];
}