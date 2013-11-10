/*
 * Copyright 2011, Blender Foundation.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Contributor:
 *		Dalai Felinto
 */

#include "COM_MapUVOperation.h"
#include "BLI_math.h"

MapUVOperation::MapUVOperation() : NodeOperation()
{
	this->addInputSocket(COM_DT_COLOR, COM_SC_NO_RESIZE);
	this->addInputSocket(COM_DT_VECTOR);
	this->addOutputSocket(COM_DT_COLOR);
	this->m_alpha = 0.0f;
	this->setComplex(true);

	this->m_inputUVProgram = NULL;
	this->m_inputColorProgram = NULL;
}

void MapUVOperation::initExecution()
{
	this->m_inputColorProgram = this->getInputSocketReader(0);
	this->m_inputUVProgram = this->getInputSocketReader(1);
}

void MapUVOperation::executePixelSampled(float output[4], float x, float y, PixelSampler sampler)
{
	float xy[2] = { x, y };
	float uv[2], deriv[2][2], alpha;

	if (!pixelTransform(xy, uv, deriv, alpha) || alpha == 0.0f) {
		zero_v4(output);
		return;
	}

	/* EWA filtering */
	this->m_inputColorProgram->readFiltered(output, uv[0], uv[1], deriv[0], deriv[1], COM_PS_NEAREST);
	
	/* UV to alpha threshold */
	const float threshold = this->m_alpha * 0.05f;
	/* XXX what is this supposed to do?!? */
	/*float alpha = 1.0f - threshold * (dx + dy);*/
	float du = len_v2(deriv[0]);
	float dv = len_v2(deriv[1]);
	float factor = 1.0f - threshold * (du + dv);
	if (factor < 0.f) alpha = 0.f;
	else alpha *= factor;

	/* "premul" */
	if (alpha < 1.0f) {
		mul_v4_fl(output, alpha);
	}
}

bool MapUVOperation::pixelTransform(const float co[2], float r_co[2], float r_deriv[2][2], float &r_alpha)
{
	float width = m_inputColorProgram->getWidth();
	float height = m_inputColorProgram->getHeight();
	float uv[4];

	m_inputUVProgram->readSampler(uv, co[0], co[1], COM_PS_NEAREST);
	r_co[0] = uv[0] * width;
	r_co[1] = uv[1] * height;
	r_alpha = uv[2];

	/* XXX currently there is no way to get real derivatives from the UV map input.
	 * Instead use a simple 1st order estimate ...
	 */
	const float epsilon[2] = { width != 0.0f ? 1.0f/width : 0.0f, height != 0.0f ? 1.0f/height : 0.0f };

	m_inputUVProgram->readSampler(uv, co[0] + epsilon[0], co[1], COM_PS_NEAREST);
	r_deriv[0][0] = uv[0];
	r_deriv[1][0] = uv[1];
	m_inputUVProgram->readSampler(uv, co[0] - epsilon[0], co[1], COM_PS_NEAREST);
	r_deriv[0][0] = 0.5f*(r_deriv[0][0] - uv[0]) * width;
	r_deriv[1][0] = 0.5f*(r_deriv[1][0] - uv[1]) * width;

	m_inputUVProgram->readSampler(uv, co[0], co[1] + epsilon[1], COM_PS_NEAREST);
	r_deriv[0][1] = uv[0];
	r_deriv[1][1] = uv[1];
	m_inputUVProgram->readSampler(uv, co[0], co[1] - epsilon[1], COM_PS_NEAREST);
	r_deriv[0][1] = 0.5f*(r_deriv[0][1] - uv[0]) * height;
	r_deriv[1][1] = 0.5f*(r_deriv[1][1] - uv[1]) * height;

	return true;
}

void MapUVOperation::deinitExecution()
{
	this->m_inputUVProgram = NULL;
	this->m_inputColorProgram = NULL;
}

bool MapUVOperation::determineDependingAreaOfInterest(rcti *input, ReadBufferOperation *readOperation, rcti *output)
{
	rcti colorInput;
	rcti uvInput;
	NodeOperation *operation = NULL;

	/* the uv buffer only needs a 3x3 buffer. The image needs whole buffer */

	operation = getInputOperation(0);
	colorInput.xmax = operation->getWidth();
	colorInput.xmin = 0;
	colorInput.ymax = operation->getHeight();
	colorInput.ymin = 0;
	if (operation->determineDependingAreaOfInterest(&colorInput, readOperation, output)) {
		return true;
	}

	operation = getInputOperation(1);
	uvInput.xmax = input->xmax + 1;
	uvInput.xmin = input->xmin - 1;
	uvInput.ymax = input->ymax + 1;
	uvInput.ymin = input->ymin - 1;
	if (operation->determineDependingAreaOfInterest(&uvInput, readOperation, output)) {
		return true;
	}

	return false;
}

