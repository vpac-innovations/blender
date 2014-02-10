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
 *		Jeroen Bakker 
 *		Monique Dewanchand
 */

#include "COM_BlurNode.h"
#include "DNA_node_types.h"
#include "COM_GaussianXBlurOperation.h"
#include "COM_GaussianYBlurOperation.h"
#include "COM_GaussianAlphaXBlurOperation.h"
#include "COM_GaussianAlphaYBlurOperation.h"
#include "COM_ExecutionSystem.h"
#include "COM_GaussianBokehBlurOperation.h"
#include "COM_FastGaussianBlurOperation.h"
#include "COM_MathBaseOperation.h"
#include "COM_SetValueOperation.h"
#include "COM_GammaCorrectOperation.h"

BlurNode::BlurNode(bNode *editorNode) : Node(editorNode)
{
	/* pass */
}

void BlurNode::convertToOperations(NodeCompiler *compiler, const CompositorContext *context) const
{
	bNode *editorNode = this->getbNode();
	NodeBlurData *data = (NodeBlurData *)editorNode->storage;
	InputSocket *inputSizeSocket = this->getInputSocket(1);
	bool connectedSizeSocket = inputSizeSocket->isConnected();

	const float size = this->getInputSocket(1)->getEditorValueFloat();
	
	CompositorQuality quality = context->getQuality();
	NodeOperation *input_operation = NULL, *output_operation = NULL;

	if (data->filtertype == R_FILTER_FAST_GAUSS) {
		FastGaussianBlurOperation *operationfgb = new FastGaussianBlurOperation();
		operationfgb->setData(data);
		operationfgb->setChunksize(context->getChunksize());
		operationfgb->setbNode(editorNode);
		
		compiler->addOperation(operationfgb);
		compiler->mapInputSocket(getInputSocket(1), operationfgb->getInputSocket(1));

		input_operation = operationfgb;
		output_operation = operationfgb;
	}
	else if (editorNode->custom1 & CMP_NODEFLAG_BLUR_VARIABLE_SIZE) {
		MathAddOperation *clamp = new MathAddOperation();
		SetValueOperation *zero = new SetValueOperation();
		zero->setValue(0.0f);
		clamp->setUseClamp(true);
		
		compiler->addOperation(clamp);
		compiler->addOperation(zero);
		compiler->mapInputSocket(getInputSocket(1), clamp->getInputSocket(0));
		compiler->addConnection(zero->getOutputSocket(), clamp->getInputSocket(1));
		
		GaussianAlphaXBlurOperation *operationx = new GaussianAlphaXBlurOperation();
		operationx->setData(data);
		operationx->setbNode(editorNode);
		operationx->setQuality(quality);
		operationx->setSize(1.0f);
		operationx->setFalloff(PROP_SMOOTH);
		operationx->setSubtract(false);
		
		compiler->addOperation(operationx);
		compiler->addConnection(clamp->getOutputSocket(), operationx->getInputSocket(0));
		
		GaussianAlphaYBlurOperation *operationy = new GaussianAlphaYBlurOperation();
		operationy->setData(data);
		operationy->setbNode(editorNode);
		operationy->setQuality(quality);
		operationy->setSize(1.0f);
		operationy->setFalloff(PROP_SMOOTH);
		operationy->setSubtract(false);
		
		compiler->addOperation(operationy);
		compiler->addConnection(operationx->getOutputSocket(), operationy->getInputSocket(0));
		
		GaussianBlurReferenceOperation *operation = new GaussianBlurReferenceOperation();
		operation->setData(data);
		operation->setbNode(editorNode);
		operation->setQuality(quality);
		
		compiler->addOperation(operation);
		compiler->addConnection(operationy->getOutputSocket(), operation->getInputSocket(1));
		
		output_operation = operation;
		input_operation = operation;
	}
	else if (!data->bokeh) {
		GaussianXBlurOperation *operationx = new GaussianXBlurOperation();
		operationx->setData(data);
		operationx->setbNode(editorNode);
		operationx->setQuality(quality);
		
		compiler->addOperation(operationx);
		compiler->mapInputSocket(getInputSocket(1), operationx->getInputSocket(1));
		
		GaussianYBlurOperation *operationy = new GaussianYBlurOperation();
		operationy->setData(data);
		operationy->setbNode(editorNode);
		operationy->setQuality(quality);

		compiler->addOperation(operationy);
		compiler->mapInputSocket(getInputSocket(1), operationy->getInputSocket(1));
		compiler->addConnection(operationx->getOutputSocket(), operationy->getInputSocket(0));

		if (!connectedSizeSocket) {
			operationx->setSize(size);
			operationy->setSize(size);
		}

		input_operation = operationx;
		output_operation = operationy;
	}
	else {
		GaussianBokehBlurOperation *operation = new GaussianBokehBlurOperation();
		operation->setData(data);
		operation->setbNode(editorNode);
		operation->setQuality(quality);
		
		compiler->addOperation(operation);
		compiler->mapInputSocket(getInputSocket(1), operation->getInputSocket(1));

		if (!connectedSizeSocket) {
			operation->setSize(size);
		}

		input_operation = operation;
		output_operation = operation;
	}

	if (data->gamma) {
		GammaCorrectOperation *correct = new GammaCorrectOperation();
		GammaUncorrectOperation *inverse = new GammaUncorrectOperation();
		graph->addOperation(correct);
		graph->addOperation(inverse);
		
		compiler->mapInputSocket(getInputSocket(0), correct->getInputSocket(0));
		compiler->addConnection(correct->getOutputSocket(), input_operation->getInputSocket(0));
		compiler->addConnection(output_operation->getOutputSocket(), inverse->getInputSocket(0));
		compiler->mapOutputSocket(getOutputSocket(), inverse->getOutputSocket());
		
		addPreviewOperation(graph, context, inverse->getOutputSocket());
	}
	else {
		compiler->mapInputSocket(getInputSocket(0), input_operation->getInputSocket(0));
		compiler->mapOutputSocket(getOutputSocket(), output_operation->getOutputSocket());
		
		addPreviewOperation(graph, context, output_operation->getOutputSocket());
	}
}
